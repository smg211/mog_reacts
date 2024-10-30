function decode_reach = DecodeReach(cfg, spike, kinematics)
% Given spiking and kinematics data, decode the reach direction from 
% spiking activity at various time lags
%
% Use as:
%       decode_reach = DecodeReach(cfg, spike, kinematics)
%
% 'spike' is an N_unit element structure containing the following fields:
%       ts: 1 x N_spike vector of times when this unit fired a spike
%
% 'kinematics' is an N_trials element structure containing the following
% fields:
%       time: 1 x N vector of times associated with each wrist velocity
%       vel_mag: 1 x N vector of the velocity vector magnitude (i.e., speed)
%       vel_theta: 1 x N vector of the velocity vector direction
%
% Sandon Griffin (2024)

% Get General Saving & Plotting Settings
doplot          = ft_getopt(cfg, 'doplot', false);

% Get function-specific settings
t_shift           = ft_getopt(cfg, 't_shift', 0:0.01:0.2); % sec
targ_sel          = ft_getopt(cfg, 'targ_sel', 'all');
docv_lags         = ft_getopt(cfg, 'docv_lags', true);
kfolds            = ft_getopt(cfg, 'kfolds', 10);
% min_samp_perclass = ft_getopt(cfg, 'min_samp_perclass', []);

n_shift = length(t_shift);

cfg = updateCFG(cfg, who);

%% LOAD DATA
n_unit = length(spike);
n_targ = size(kinematics.trial(1).t_targon, 1);
if strcmp(targ_sel, 'all')
  targ_sel = 1:n_targ;
end
cfg.targ_sel = targ_sel;

i_trl = 1:length(kinematics.trial);
n_trl = length(kinematics.trial);

% INITIALIZE STRUCT FOR OUTPUT
decode_reach = [];
decode_reach.cfg = cfg;

%% DETERMINE THE SPIKE RATE OF EACH NEURON DURING EACH REACH DIRECTION
% first get a raster of all spiking activity
t_bin_rast = 0.001;
[n_spk_all, t_all] = make_raster_stride(spike, t_bin_rast, t_bin_rast);

nbin_shift = round(t_shift./t_bin_rast);

cfg_toi = cfg;
cfg_toi.kinematics = kinematics;
cfg_toi.trial_sel = i_trl;
cfg_toi.targ_sel = 2:5; % only use the time from when targets 2 to 5 are on the screen, since they are optimally arranged to have one reach direction in each quadrant
cfg_toi.topsec = 2.5; % use the bins corresponding to the 2.5 seconds where speed is the highest in each direction
cfg_toi.binsize_theta = pi/4; % divide reach direction into 8 quadrants
cfg_toi.stride_theta = pi/4; % make the quadrants non-overlapping


[is_toi_class, dur_class, class_info] = get_is_class(cfg_toi, t_all);
nclass = size(is_toi_class, 1);

fxmatrix = nan(n_shift, n_unit, nclass);
for t = 1:length(t_shift)
  for class = 1:nclass
    is_toi_shift = circshift(is_toi_class(class, :), -nbin_shift(t));
    
    for n = 1:n_unit
      n_spk_targ_n = sum(n_spk_all(n, is_toi_shift));
      fxmatrix(t, n, class) = n_spk_targ_n/dur_class(class);
    end
  end
end

fxmatrix(isnan(fxmatrix)) = 0;
% assert(~any(isnan(fxmatrix(:))));

% store output
decode_reach.fxmatrix = fxmatrix;
decode_reach.class_info = class_info;
decode_reach.dur_class = dur_class;

%% BIN THE SPIKES FOR DECODING AND GET THETA REACH DURING EACH BIN
binsize_decode = 0.1; % sec
stride_decode = 0.02; % sec
[~, t_decode, t_binedges] = make_raster_stride(spike, binsize_decode, stride_decode);

nbin_decode = length(t_decode);
cfg_toi.doplot = false;
[is_toi_decode_all, ~, ~] = get_is_class(cfg_toi, t_decode);
is_toi_decode_any = any(is_toi_decode_all, 1);

% get the average reach direction during each bin and then assign that
% bin to the theta bin with the closest center point to that average
t_k = [kinematics.trial.time];
theta = [kinematics.trial.vel_theta];
theta_bincents = class_info.cond_vals;

theta_true_raw = nan(1, nbin_decode);
class_true_ix = nan(1, nbin_decode);
for t = 1:size(t_binedges, 1)
  if any(t_k > t_binedges(t, 1) & t_k <= t_binedges(t, 2))
    i_t_k = find(t_k > t_binedges(t, 1) & t_k <= t_binedges(t, 2));
    theta_true_raw(t) = circ_mean(theta(i_t_k), [], 2);
    [~, class_true_ix(t)] = min(abs(theta_bincents - theta_true_raw(t)));
  end
end

% store output
decode_reach.is_toi_decode_any = is_toi_decode_any;
decode_reach.class_true = class_true_ix;
decode_reach.class_true_raw = theta_true_raw;

%% PERFORM DECODING BEFORE CROSS-VALIDATION
% AT EACH DIFFERENT LAG, CALCULATE P(class)
p_class_pred = nan(n_shift, nclass, sum(is_toi_decode_any));
class_pred_ix = nan(n_shift, sum(is_toi_decode_any));
for t = 1:n_shift
  % shift the spikes
  spike_shift = spike;
  for u = 1:length(spike)
    spike_shift(u).ts = spike(u).ts + t_shift(t);
  end
  [n_spk_decode, ~, ~] = make_raster_stride(spike_shift, binsize_decode, stride_decode);
  
  n_spk_decode_sel = n_spk_decode(:, is_toi_decode_any);
  
  [p_class_pred(t, :, :), class_pred_ix(t, :)] = do_bayes_decode(squeeze(fxmatrix(t, :, :)), ...
    n_spk_decode_sel, binsize_decode);
end

%% CROSS-VALIDATE DECODER ACCURACY AT EACH LAG USING K-FOLD CROSS VALIDATION
if docv_lags
  % randomly split the trials into k groups
  assert(rem(n_trl, kfolds) == 0);
  n_trlpergrp = n_trl/kfolds;
  i_grplims = round(linspace(0, n_trl, kfolds+1));
  trlorder = randperm(n_trl);
  
  i_trl_grp = nan(kfolds, n_trlpergrp);
  for k = 1:kfolds
    i_trl_sel = trlorder((i_grplims(k)+1):i_grplims(k+1));
    i_trl_grp(k, :) = i_trl(i_trl_sel);
  end
  
  % for each fold, get fxmatrix only using the train data, and then predict
  % direction on the held out trials
  doinclude_class_cv = true(nclass, kfolds);
  p_class_pred_cv = cell(1, kfolds);
  
  class_pred_ix_cv = cell(1, kfolds);
  class_pred_test_all_cv = cell(1, n_shift);
  p_class_pred_all_cv = cell(1, n_shift);
  class_true_test_all = [];
  theta_true_test_raw_all = [];
  for t = 1:n_shift
    class_pred_test_all_cv{t} = [];
    p_class_pred_all_cv{t} = [];
  end
  tic
  for k = 1:kfolds
    i_k = [1:k-1 k+1:kfolds];
    i_trl_k = i_trl_grp(i_k, :);
    
    % get the indices of times in the all spike raster associated with the
    % training data
    cfg_toi_train = cfg_toi;
    cfg_toi_train.trial_sel = i_trl_k(:);
    cfg_toi_train.doplot = false;
    cfg_toi_train.speed_thresh = class_info.speed_thresh;
    [is_toi_train, dur_class] = get_is_class(cfg_toi_train, t_all);
    
    
    % get the indices of times in the decode raster associated with the held
    % out trials
    cfg_toi_test = cfg_toi;
    cfg_toi_test.trial_sel = i_trl_grp(k, :);
    cfg_toi_test.doplot = false;
    cfg_toi_test.speed_thresh = class_info.speed_thresh;
    is_toi_decode_test = get_is_class(cfg_toi_test, t_decode);
    is_toi_test_any = any(is_toi_decode_test, 1);
    
    nbin_decode_test = sum(is_toi_test_any);
    
    class_true_test = class_true_ix(is_toi_test_any);
    class_true_test_all = [class_true_test_all class_true_test];
    
    theta_true_test_raw = theta_true_raw(is_toi_test_any);
    theta_true_test_raw_all = [theta_true_test_raw_all theta_true_test_raw];
    
    % AT EACH DIFFERENT LAG, CALCULATE P(class)
    p_class_pred_cv{k} = nan(n_shift, nclass, nbin_decode_test);
    class_pred_ix_cv{k} = nan(n_shift, nbin_decode_test);
    fxmatrix_cv = nan(n_shift, n_unit, nclass);
    for t = 1:n_shift
      % get the tuning profiles for each classition
      for class = 1:nclass
        is_toi_shift = circshift(is_toi_train(class, :), -nbin_shift(t));
        
        for n = 1:n_unit
          n_spk_targ_n = sum(n_spk_all(n, is_toi_shift));
          fxmatrix_cv(t, n, class) = n_spk_targ_n/dur_class(class);
        end
      end
      
      fxmatrix_cv(isnan(fxmatrix_cv)) = 0;
      
      % shift the spikes
      spike_shift = spike;
      for u = 1:length(spike)
        spike_shift(u).ts = spike(u).ts + t_shift(t);
      end
      [n_spk_decode_shift, ~, ~] = make_raster_stride(spike_shift, ...
        binsize_decode, stride_decode);
      
      n_spk_decode_test = n_spk_decode_shift(:, is_toi_test_any);
      
      [p_class_pred_cv{k}(t, :, :), class_pred_ix_cv{k}(t, :)] = do_bayes_decode(squeeze(fxmatrix_cv(t, :, :)), ...
        n_spk_decode_test, binsize_decode);
      
      class_pred_test_all_cv{t} = [class_pred_test_all_cv{t} class_pred_ix_cv{k}(t, :)];
      p_class_pred_all_cv{t} = [p_class_pred_all_cv{t} squeeze(p_class_pred_cv{k}(t, :, :))];
    end
  end
  
  
  % store output
  decode_reach.doinclude_class_cv = doinclude_class_cv;
  decode_reach.p_class_pred_cv = p_class_pred_cv;
  decode_reach.class_pred_ix_cv = class_pred_ix_cv;
  
  %% CALCULATE DECODING ACCURACY FOR EACH TIME LAG
  conf_mat_cv = nan(n_shift, nclass, nclass);
  p_conf_mat_cv = nan(n_shift, nclass, nclass);
  rho_lag_cv = nan(1, n_shift);
  sse_lag_cv = nan(1, n_shift);
  error_lag_cv = nan(n_shift, length(class_true_test_all));
  for t = 1:n_shift
    % get the confusion matrix and percent correct, and error for each time lag
    for i = 1:nclass
      true_is_i = class_true_test_all == i;
      for j = 1:nclass
        conf_mat_cv(t, i, j) = sum(class_pred_test_all_cv{t}(true_is_i) == j);
      end
      p_conf_mat_cv(t, i, :) = conf_mat_cv(t, i, :)./sum(true_is_i);
    end
    
    % get the error
    theta_true_toi = theta_bincents(class_true_test_all);
    theta_pred_toi = theta_bincents(class_pred_test_all_cv{t});
    
    rho_lag_cv(t) = circ_corrcc(theta_true_toi, theta_pred_toi);
    error_lag_cv(t, :) = angdiff(theta_true_toi, theta_pred_toi);
    sse_lag_cv(t) = sum(error_lag_cv(t, :).^2);
  end
  
  [~, i_lag_sel] = min(sse_lag_cv);
  t_shift_sel = t_shift(i_lag_sel);
  
  % PLOT THE NORMALIZED ERROR AND RHO
  if doplot
    fig = figure; hold on;
    yyaxis left
    plot(t_shift, rho_lag_cv.^2, '-b', 'LineWidth', 3);
    
    yyaxis right
    plot(t_shift, sse_lag_cv./max(sse_lag_cv), '-r', 'LineWidth', 3);
    
    a = gca;
    ylims = a.YAxis(2).Limits;
    
    plot(t_shift_sel([1 1]), ylims, '--k');
    
    a.YAxis(2).Limits = ylims;
    
    a.XLabel.String = 'Lags (>0 = neural leading behavior';
    a.YAxis(2).Label.String = 'Sum of Squared Errors';
    a.YAxis(1).Label.String = 'R-squared';
    a.FontSize = 18;
  end
  
  % store output
  decode_reach.conf_mat_cv = conf_mat_cv;
  decode_reach.p_conf_mat_cv = p_conf_mat_cv;
  decode_reach.rho_lag_cv = rho_lag_cv;
  decode_reach.sse_lag_cv = sse_lag_cv;
end
