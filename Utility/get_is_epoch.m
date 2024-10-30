function [is_toi, toi_dur] = get_is_epoch(cfg, t, toi)
%% Get a 1 x T logical array of T times in t that occur during the toi
% note a lot of the code is copied from the code I used for the flint
% dataset
%
% Use as: 
%       [is_toi, toi_dur] = get_is_epoch(cfg, t, toi)
%
% Sandon Griffin (2024)

% Get inputs
kinematics          = ft_getopt(cfg, 'kinematics', []);
speed_thresh        = ft_getopt(cfg, 'speed_thresh', []);
theta_binedges      = ft_getopt(cfg, 'theta_binedges', []);

%% MAIN
toi_dur = [];


if isempty(t)
  is_toi = [];
elseif strcmp(toi, 'all')
  is_toi = logical(ones(1, length(t)));
else
  % only take the bins that correspond to times that we are training on
  if any(strfind(toi, '_and_')) 
    toi_operator = 'and';
    i_op = strfind(toi, '_and_');
    i_last = 1;
    for i = 1:length(i_op) + 1
      if i == length(i_op) + 1
        toi_list{i} = toi(i_op(i-1)+5:end);
      elseif i == 1
        toi_list{i} = toi(1:i_op(i)-1);
      else
        toi_list{i} = toi(i_op(i-1)+5:i_op(i)-1);
      end
    end
  elseif any(strfind(toi, '_or_'))
    toi_operator = 'or';
    i_op = strfind(toi, '_or_');
    i_last = 1;
    for i = 1:length(i_op) + 1
      if i == length(i_op) + 1
        toi_list{i} = toi(i_op(i-1)+4:end);
      elseif i == 1
        toi_list{i} = toi(1:i_op(i)-1);
      else
        toi_list{i} = toi(i_op(i-1)+4:i_op(i)-1);
      end
    end
  else
    toi_list{1} = toi;
    toi_operator = [];
  end

  is_toi = [];
  for i_toi = 1:length(toi_list)
    % DETERMINE IF THIS TOI SHOULD ONLY BE APPLIED TO SOME TARGET INDEXES
    [toi_list{i_toi}, targ_ix] = get_targ_ix(toi_list{i_toi});

    % TIMES DURING BREAK PERIODS
    if strfind(toi_list{i_toi}, 'break') == 1
      % determine if there is a reference to a specific break
      if length(toi_list{i_toi}) >= 6 && ~isempty(str2num(toi_list{i_toi}(6:end))) && ...
          isnumeric(str2num(toi_list{i_toi}(6:end)))
        break_ix = str2num(toi_list{i_toi}(6:end));
        is_toi(i_toi, :) = any(t >= [kinematics.break(break_ix).t_start]' & t <= [kinematics.break(break_ix).t_end]', 1);
        
        toi_dur = nansum([kinematics.break(break_ix).t_end] - [kinematics.break(break_ix).t_start]);
      elseif isfield(kinematics, 'break')
        is_toi(i_toi, :) = any(t >= [kinematics.break.t_start]' & t <= [kinematics.break.t_end]', 1);
        toi_dur = nansum([kinematics.break.t_end] - [kinematics.break.t_start]);
      else
        is_toi(i_toi, :) = false(1, length(t));
        toi_dur = 0;
      end
      
    % TIMES DURING BREAK MOVEMENTS
    elseif strcmp(toi_list{i_toi}, 'brkmove')
      if length([kinematics.break.t_move_start]) > 0
        is_toi(i_toi, :) = any(t >= [kinematics.break.t_move_start]' & ...
          t <= [kinematics.break.t_move_end]', 1);
        toi_dur = nansum([kinematics.break.t_move_start] - [kinematics.break.t_move_end]);
      else
        is_toi(i_toi, :) = false(1, length(t));
        toi_dur = 0;
      end
      
      
    % TIMES DURING TARGET ON THE SCREEN
    elseif strfind(toi_list{i_toi}, 'targon') == 1
      t_on = [kinematics.trial.t_targon];
      t_off = [kinematics.trial.t_targoff];
      if isempty(targ_ix)
        t_on = min(t_on, [], 1);
        t_off = max(t_off, [], 1);
      else
        t_on = t_on(targ_ix, :);
        t_off = t_off(targ_ix, :);
      end
      
      is_toi(i_toi, :) = any(t >= t_on' & t <= t_off', 1);
      toi_dur = nansum(t_off - t_on);

    % TIMES DURING TRIALS OF A SPECIFIC NUMBER
    elseif any(strfind(toi_list{i_toi}, 'trial'))
      if isfield(cfg, 'trial_sel')
        i_trial = cfg.trial_sel;
      else
        i_trial = str2num(toi_list{i_toi}(6:length(toi_list{i_toi})));
      end
      is_toi(i_toi, :) = any(t >= [kinematics.trial(i_trial).t_start]' & t <= [kinematics.trial(i_trial).t_end]', 1);
      toi_dur = nansum([kinematics.trial(i_trial).t_end] - [kinematics.trial(i_trial).t_start]);
      
    % TIMES DURING BLOCKS OF A SPECIFIC NUMBER
    elseif any(strfind(toi_list{i_toi}, 'block'))
      i_block = str2num(toi_list{i_toi}(6:length(toi_list{i_toi})));
      i_1sttrial = find([kinematics.trial.block_ix] == i_block, 1, 'first');
      i_lasttrial = find([kinematics.trial.block_ix] == i_block, 1, 'last');
      t_start = kinematics.trial(i_1sttrial).t_start;
      t_end = kinematics.trial(i_lasttrial).t_end;
      
      is_toi(i_toi, :) = any(t >= t_start & t <= t_end, 1);
      toi_dur = nansum(t_end - t_start);
      
    % TIMES DURING REACHES IN A GIVEN DIRECTION
    elseif strcmp(toi_list{i_toi}, 'reachdir')      
      t_theta_start = [];
      t_theta_end = [];
      for trl = 1:length(kinematics.trial)
        t_trl = kinematics.trial(trl).time;
        theta_trl = kinematics.trial(trl).vel_theta;
        is_theta = theta_trl > theta_binedges(1) & theta_trl <= theta_binedges(2);
        
        i_start = find(diff([0 is_theta]) == 1);
        i_end = find(diff([is_theta 0]) == -1);
        t_theta_start = [t_theta_start t_trl(i_start)];
        t_theta_end = [t_theta_end t_trl(i_end)];
      end
      
      
      is_toi(i_toi, :) = any(t > t_theta_start' & t <= t_theta_end', 1);
      toi_dur = sum(t_theta_end - t_theta_start);
      
    % TIMES DURING REACHES ABOVE A GIVEN SPEED
    elseif strcmp(toi_list{i_toi}, 'reachspeed')
      t_speed_start = [];
      t_speed_end = [];
      for trl = 1:length(kinematics.trial)
        t_trl = kinematics.trial(trl).time;
        speed_trl = kinematics.trial(trl).vel_mag;
        is_abv_speed_thresh = speed_trl > speed_thresh;
        
        i_start = find(diff([0 is_abv_speed_thresh]) == 1);
        i_end = find(diff([is_abv_speed_thresh 0]) == -1);
        t_speed_start = [t_speed_start t_trl(i_start)];
        t_speed_end = [t_speed_end t_trl(i_end)];
      end
      
      
      is_toi(i_toi, :) = any(t >= t_speed_start' & t <= t_speed_end', 1);
      toi_dur = sum(t_speed_end - t_speed_start);
    else
      error(['This function does not know how to deal with the time of interest: ' toi_list{i_toi}]);
    end
  end
  if strcmp(toi_operator, 'and') || isempty(toi_operator)
    is_toi = all(is_toi, 1);
  elseif strcmp(toi_operator, 'or')
    is_toi = any(is_toi, 1);
  end
end



%% SUBFUNCTIONS
function [field, targ_ix] = get_targ_ix(input)
field = input;
targ_ix = [];
if any(strfind(input, 'targ'))
  i_targ_str = strfind(input, 'targ');
  if length(input) >= i_targ_str+4 && ~isempty(str2num(input(i_targ_str+4))) && ...
      isnumeric(str2num(input(i_targ_str+4)))
    targ_ix = str2num(input(i_targ_str+4));
    field = [input(1:i_targ_str+3) input(i_targ_str+5:end)];
  end
end

