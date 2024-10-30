function [kinematics, spike, rr_brk_targ_gt] = SimData(cfg)
% Create simulated reach trajectories and spike data during tasks and breaks

spike2kin_lag     = ft_getopt(cfg, 'spike2kin_lag', 0);

n_trl = 40;
n_trl_perblock = 10;

%% set target locations and get angle between targets
targ = [0 -0.065; ...
  0.065 0.065; ...
  -0.065 0.065; ...
  -0.065 -0.065; ...
  0.065 -0.065];

theta_targ2targ = nan(1, 5);
for t = 1:5
  if t == 1
    theta_targ2targ(t) = atan2(targ(t, 2)-(-0.15), targ(t, 1));
  else
    theta_targ2targ(t) = atan2(targ(t, 2)-targ(t-1, 2), targ(t, 1)-targ(t-1, 1));
  end
end

%% generate wrist trajectories towards those targets
% give every movement the same velocity profile
targ_rad = 0.019;
dt_k = 0.02;

decel_factor = 0.95;
accel_factor = 1.05;
noise_factor = 0.0005;
sm_kernel = 20;

kinematics = [];
kinematics.theta_targ2targ = theta_targ2targ;

% figure;
for trl = 1:n_trl
  if trl == 1
    kinematics.trial(trl).t_start = 0;
  else
    if rem(trl, n_trl_perblock) == 1
      i_brk = floor(trl/n_trl_perblock);
      kinematics.break(i_brk).t_start = kinematics.trial(trl-1).t_end + 2;
      kinematics.break(i_brk).t_end = kinematics.trial(trl-1).t_end + 92;
      kinematics.trial(trl).t_start = kinematics.trial(trl-1).t_end + 92;
    else
      kinematics.trial(trl).t_start = kinematics.trial(trl-1).t_end + 2;
    end
  end
  
  wr_pos = [0 -0.15; 0 -0.145];
  ix = 2;
  ix_targon = nan(5, 1);
  ix_targoff = nan(5, 1);
  for t = 1:5
    ix_targon(t) = ix;
    targon = targ(t, :);
    d2targ_total = sqrt(sum((wr_pos(end, :)-targon).^2));
    while sqrt(sum((wr_pos(end, :)-targon).^2)) > targ_rad
      ix = ix + 1;
      % get the distance and angle to the target
      d2targ = sqrt(sum((wr_pos(end, :)-targon).^2));
      theta2targ = atan2(targon(2)-wr_pos(end, 2), targon(1)-wr_pos(end, 1));
      dxy = sqrt(sum(diff(wr_pos(end-1:end, :)).^2));
      d_wr_last_x = diff(wr_pos(end-1:end, 1));
      d_wr_last_y = diff(wr_pos(end-1:end, 2));
      wr_theta = atan2(wr_pos(end, 2)-wr_pos(end-1, 2), wr_pos(end, 1)-wr_pos(end-1, 1));
      
      % decelerate to a stop and then accelerate towards the target
      if abs(angdiff(theta2targ, wr_theta)) > pi/4
        if dxy > 0.002
          % wrist is still moving with some speed
          decel_accel = -1;
        else
          % wrist has slowed almost to a stop
          d2targ_total = d2targ;
          decel_accel = 1;
        end
      else
        if d2targ < 0.5*d2targ_total
          decel_accel = -1;
        else
          decel_accel = 1;
        end
      end
      
      if decel_accel == -1
        d_wr_next_x = cos(wr_theta)*dxy*decel_factor;
        d_wr_next_y = sin(wr_theta)*dxy*decel_factor;
        wr_pos(ix, 1) = wr_pos(ix-1, 1) + d_wr_next_x + noise_factor*rand-(noise_factor/2);
        wr_pos(ix, 2) = wr_pos(ix-1, 2) + d_wr_next_y + noise_factor*rand-(noise_factor/2);
      elseif decel_accel == 1
        d_wr_next_x = cos(theta2targ)*dxy*accel_factor;
        d_wr_next_y = sin(theta2targ)*dxy*accel_factor;
        wr_pos(ix, 1) = wr_pos(ix-1, 1) + d_wr_next_x + noise_factor*rand-(noise_factor/2);
        wr_pos(ix, 2) = wr_pos(ix-1, 2) + d_wr_next_y + noise_factor*rand-(noise_factor/2);
      end
    end
    ix_targoff(t) = ix;
  end
  
  kinematics.trial(trl).time = (dt_k:dt_k:(dt_k*ix)) + kinematics.trial(trl).t_start;
  kinematics.trial(trl).t_end = kinematics.trial(trl).time(end) + 0.5;
  kinematics.trial(trl).t_targon = kinematics.trial(trl).time(ix_targon)';
  kinematics.trial(trl).t_targoff = kinematics.trial(trl).time(ix_targoff)';
  kinematics.trial(trl).pos = smooth_dat(wr_pos', sm_kernel, 'moving')';
  kinematics.trial(trl).vel_mag = [sqrt(sum(diff(kinematics.trial(trl).pos', 1, 2).^2, 1)) nan];
  kinematics.trial(trl).vel_theta = [atan2(diff(kinematics.trial(trl).pos(:, 2)', 1, 2), diff(kinematics.trial(trl).pos(:, 1)', 1, 2)) nan];
  
%   subplot(1, 2, 1); hold on;
%   plot(kinematics.trial(trl).pos(:, 1), kinematics.trial(trl).pos(:, 2));
%   a = gca;
%   a.YLim = [-0.2 0.2];
%   axis equal
%   
%   subplot(1, 2, 2); hold on;
%   plot(kinematics.trial(trl).vel_mag)
end

%% generate tuning profiles for N neurons
n_dir = 100;
n_unit = 100;


% each unit should have a tuning angle, tuning width, and tuning strength
dir = linspace(-pi, pi, n_dir+1);
dir = dir(2:end);

fxmatrix = nan(n_unit, n_dir);
for u = 1:n_unit
  fx_u = zeros(1, n_dir);
  
  theta_tune = dir(u);
  tune_width = 1;
  tune_strength = 1;

  y = (cos(tune_width*dir)+1)/2;
  y = (1-tune_strength)+(y*tune_strength);
  fx_u(:) = (1-tune_strength);
  iy = (n_dir/2)+(-round((n_dir/(tune_width*2))):round((n_dir/(tune_width*2))));
  iy(iy < 1) = [];
  fx_u(iy) = y(iy);
  fx_u = circshift(fx_u, round(n_dir*theta_tune/(2*pi)));
  fxmatrix(u, :) = fx_u;
end

%% generate spiking data according to the tuning profiles and velocity
dt_spike = 0.001;
max_vel_mag = max([kinematics.trial.vel_mag]);

for u = 1:n_unit
  spike(u).ts = [];
end

for trl = 1:n_trl
  % upsample the velocity data
  t_new = linspace(kinematics.trial(trl).time(1), kinematics.trial(trl).time(end), round(length(kinematics.trial(trl).time)*dt_k/dt_spike));
  tsin = timeseries(kinematics.trial(trl).vel_mag, kinematics.trial(trl).time);
  tsout = resample(tsin, t_new);
  mag_us = squeeze(tsout.Data);
  mag_us_norm = mag_us./(0.5*max_vel_mag);
  mag_us_norm(mag_us_norm > 1) = 1;
  
  tsin = timeseries(kinematics.trial(trl).vel_theta, kinematics.trial(trl).time);
  tsout = resample(tsin, t_new);
  theta_us = squeeze(tsout.Data);
  
  [~, i_theta] = min(abs(dir - theta_us), [], 2);
  
  p_spike = fxmatrix(:, i_theta).*mag_us_norm'.*poissrnd(1, [n_unit, length(t_new)])./5;
  rast = p_spike > 0.3;
  for u = 1:n_unit
    t_spk = t_new(find(rast(u, :))) - spike2kin_lag;
    spike(u).ts = [spike(u).ts t_spk];
  end
end

%% generate spiking data during breaks with short bouts of reactivations
rr_brk_gt = 1+sort(rand(1, length(kinematics.break))-0.5, 'descend');
rr_brk_targ_gt = nan(5, length(kinematics.break));

fxmatrix(:, 101) = 1;

[~, i_dir_reach] = min(abs(dir - theta_targ2targ'), [], 2);

for k = 1:length(kinematics.break)
  kinematics.break(k).t_move_start = [];
  kinematics.break(k).t_move_end = [];
  
  rr_targ = 1+rand(1, 5)-0.5;
  t_brk = kinematics.break(k).t_start:0.02:kinematics.break(k).t_end;
  i_theta = ones(1, length(t_brk))+length(dir);
  mag_react = 0.35*ones(1, length(t_brk));
  for t = 1:5
    rr_brk_targ_gt(t, k) = rr_brk_gt(k)*rr_targ(t);
    nreact_kt = round(rr_brk_targ_gt(t, k)*90);
    
    i_react = randperm(length(t_brk), nreact_kt);
    i_theta(i_react) = i_dir_reach(t);
    mag_react(i_react) = 0.9;
  end
  
  t_new = linspace(kinematics.break(k).t_start, kinematics.break(k).t_end, round(length(t_brk)*dt_k/dt_spike));
  tsin = timeseries(mag_react, t_brk);
  tsout = resample(tsin, t_new);
  mag_react_us = squeeze(tsout.Data);
  
  tsin = timeseries(i_theta, t_brk);
  tsout = resample(tsin, t_new);
  i_theta_us = round(squeeze(tsout.Data));
  
  p_spike = fxmatrix(:, i_theta_us).*mag_react_us'.*poissrnd(1, [n_unit, length(t_new)])./5;
  rast = p_spike > 0.3;
  for u = 1:n_unit
    t_spk = t_new(find(rast(u, :)));
    spike(u).ts = [spike(u).ts t_spk];
  end
end

for u = 1:n_unit
  spike(u).ts = sort(spike(u).ts);
end

% store the ground truth reactivation rates
kinematics.rr_brk_targ_gt = rr_brk_targ_gt;








