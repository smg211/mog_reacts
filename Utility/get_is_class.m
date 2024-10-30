function [is_toi_cond, dur_cond, class_info] = get_is_class(cfg, t_all)
binsize_theta       = ft_getopt(cfg, 'binsize_theta', pi/4);
stride_theta        = ft_getopt(cfg, 'stride_theta', pi/4);
speed_thresh        = ft_getopt(cfg, 'speed_thresh', []);
topsec              = ft_getopt(cfg, 'topsec', 2.5);
targ_sel            = ft_getopt(cfg, 'targ_sel', []);

if isfield(cfg, 'speed_thresh')
  get_new_speed_thresh = false;
else
  get_new_speed_thresh = true;
end

t_bin_rast = mode(diff(t_all));

class_info = [];
class_info.cond_str = 'reach directions';

% get velocity magnitude and direction at all times
t_k = [cfg.kinematics.trial.time];
speed = [cfg.kinematics.trial.vel_mag];
theta = [cfg.kinematics.trial.vel_theta];

% get logical matrix of selected targets
if ~isempty(targ_sel)
  is_targon_spk = [];
  is_targon_k = [];
  for t = 1:length(targ_sel)
    is_targon_spk(t, :) = get_is_epoch(cfg, t_all, ['targ' num2str(targ_sel(t)) 'on']);
    is_targon_k(t, :) = get_is_epoch(cfg, t_k, ['targ' num2str(targ_sel(t)) 'on']);
  end
  is_targon_spk = any(is_targon_spk, 1);
  is_targon_k = any(is_targon_k, 1);
else
  is_targon_spk = get_is_epoch(cfg, t_all, 'targon');
  is_targon_k = get_is_epoch(cfg, t_k, 'targon');
end

% get logical matrix of selected trials
if isfield(cfg, 'trial_sel')
  is_trial_sel = get_is_epoch(cfg, t_all, 'trial');
  is_trial_sel_k = get_is_epoch(cfg, t_k, 'trial');
else
  is_trial_sel = true(1, length(t_all));
  is_trial_sel_k = true(1, length(t_k));
end

% get binedges and bincenters for each direction
n_thetabins = round(2*pi/(stride_theta));
theta_binedges(1, :) = [0 binsize_theta];
for b = 2:n_thetabins
  theta_binedges(b, :) = theta_binedges(b-1, :)+stride_theta;
end
theta_binedges = theta_binedges+binsize_theta/2;
theta_binedges = rem(theta_binedges, 2*pi);
theta_binedges(theta_binedges > pi) = theta_binedges(theta_binedges > pi)-2*pi;
theta_bincents = circ_mean(theta_binedges, [], 2);

% get logical matrix of selected directions where speed exceeds the given
% threshold
hc_binedges = linspace(-pi, pi, 32);
hc_bincents = get_bincents(hc_binedges);
hc_theta = nan(n_thetabins, length(hc_bincents));
is_toi_cond = false(n_thetabins, length(t_all));
is_toi_kintime = false(n_thetabins, length(t_k));
dur_cond = nan(1, n_thetabins);
for b = 1:n_thetabins
  if sign(theta_binedges(b, 1)) ~= sign(theta_binedges(b, 2)) && ~any(theta_binedges(b, :) == 0)
    i_neg = find(theta_binedges(b, :) < 0);
    i_pos = find(theta_binedges(b, :) > 0);
    [~, is_zero_or_pi] = min(abs([0 pi] - nanmean(abs(theta_binedges(b, :)))));
    if is_zero_or_pi == 1
      cfg.theta_binedges = [theta_binedges(b, i_neg) 0];
      is_theta_neg = get_is_epoch(cfg, t_all, 'reachdir');
      is_theta_neg_k = get_is_epoch(cfg, t_k, 'reachdir');
      
      cfg.theta_binedges = [0 theta_binedges(b, i_pos)];
      is_theta_pos = get_is_epoch(cfg, t_all, 'reachdir');
      is_theta_pos_k = get_is_epoch(cfg, t_k, 'reachdir');
    elseif is_zero_or_pi == 2
      cfg.theta_binedges = [-pi theta_binedges(b, i_neg)];
      is_theta_neg = get_is_epoch(cfg, t_all, 'reachdir');
      is_theta_neg_k = get_is_epoch(cfg, t_k, 'reachdir');
      
      cfg.theta_binedges = [theta_binedges(b, i_pos) pi];
      is_theta_pos = get_is_epoch(cfg, t_all, 'reachdir');
      is_theta_pos_k = get_is_epoch(cfg, t_k, 'reachdir');
    end
    
    is_theta = is_theta_neg | is_theta_pos;
    is_theta_k = is_theta_neg_k | is_theta_pos_k;
  else
    cfg.theta_binedges = theta_binedges(b, :);
    is_theta = get_is_epoch(cfg, t_all, 'reachdir');
    is_theta_k = get_is_epoch(cfg, t_k, 'reachdir');
  end
  
  if get_new_speed_thresh
    speed_targon_theta = sort(speed(is_targon_k & is_theta_k & is_trial_sel_k), 'descend');
    nbin_sec = ceil(topsec/mode(diff(t_k)));
    nbin_sec = min([nbin_sec length(speed_targon_theta)]);
    if nbin_sec == 0
      speed_thresh(b) = nan;
    else
      speed_thresh(b) = speed_targon_theta(nbin_sec);
      
      cfg.speed_thresh = speed_thresh(b);
      is_abv_speed_thresh = get_is_epoch(cfg, t_all, 'reachspeed');
      is_abv_speed_thresh_k = get_is_epoch(cfg, t_k, 'reachspeed');
    end
  else
    cfg.speed_thresh = speed_thresh(b);
    is_abv_speed_thresh = get_is_epoch(cfg, t_all, 'reachspeed');
    is_abv_speed_thresh_k = get_is_epoch(cfg, t_k, 'reachspeed');
  end
  
  is_toi_cond(b, :) = is_targon_spk & is_abv_speed_thresh & is_theta & is_trial_sel;
  dur_cond(b) = sum(is_toi_cond(b, :))*t_bin_rast;
  is_toi_kintime(b, :) = is_targon_k & is_abv_speed_thresh_k & is_theta_k & is_trial_sel_k;
  
  hc_theta(b, :) = histcounts(theta(is_toi_kintime(b, :)), hc_binedges);
end

class_info.cond_val_str = '\Theta Bincenter';
class_info.cond_vals = theta_bincents;

theta_bincents_pos = theta_bincents;
theta_bincents_pos(round(theta_bincents_pos, 2) <= 0) = ...
  theta_bincents_pos(round(theta_bincents_pos, 2) <= 0) + 2*pi;
class_info.cond_vals_pos = theta_bincents_pos;
class_info.theta_binedges = theta_binedges;
class_info.speed_thresh = speed_thresh;
class_info.is_toi_kintime = is_toi_kintime;
