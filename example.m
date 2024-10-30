[path,name,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(path));

%% CREATE SIMULATED DATA
spike2kin_lag = 0.1;

cfg = [];
cfg.spike2kin_lag = 0.1;

[kinematics, spike] = SimData(cfg);

%% DETECT REACTIVATIONS DURING BREAKS
cfg = [];
cfg.dorandperm          = false;
cfg.nperms              = 1000;

react = PCAReact(cfg, spike, kinematics);

% Plot the detected break reactivation rates vs. the ground truth reach
% reactivation rates
figure;
plot(nanmean(kinematics.rr_brk_targ_gt, 1), nanmean(react.rate_pc_brk, 1), '.k', 'MarkerSize', 40);
drawnow;
a = gca;
a.XLabel.String = 'Break React Rate (Ground Truth)';
a.YLabel.String = 'Break React Rate (PCA Detected)';

%% DECODE REACHES FROM SPIKING ACTIVITY
cfg = [];
cfg.doplot              = 1;
cfg.t_shift             = spike2kin_lag + (-0.05:0.01:0.05);
cfg.kfolds              = 10;
cfg.docv_lags           = true;

decode_reach = DecodeReach(cfg, spike, kinematics);

%% BAYESIAN DECODING OF REACTIVATIONS
cfg = [];
cfg.lag_sel             = spike2kin_lag;
cfg.nperms              = 1000;
cfg.dorandperm          = false;

decode_react = DecodeReact(cfg, spike, kinematics, decode_reach);


% Plot the detected reach reactivation rates for each break vs. the ground truth reach
% reactivation rates for each break
[~, i_min] = min(abs(decode_reach.class_info.cond_vals' - kinematics.theta_targ2targ'), [], 2);
rr_dir_brk = decode_react.rate_dir_brk(i_min, :);
figure; hold on;
cmap = cbrewer('qual', 'Set1', 5);
for t = 2:5
  plot(kinematics.rr_brk_targ_gt(t, :), rr_dir_brk(t, :), '.', 'Color', cmap(t, :), 'MarkerSize', 20);
end
a = gca;
a.XLabel.String = 'TargetBreak React Rate (Ground Truth)';
a.YLabel.String = 'TargetBreak React Rate (Detected)';
a.FontSize = 18;
leg = legend;
leg.String = {'Reach 2', 'Reach 3', 'Reach 4', 'Reach 5'};
leg.Location = 'EastOutside';













