function react = PCAReact(cfg, spike, kinematics)
% Calculate reactivation strengths and detect reactivation events using a 
% PCA-based method (first described in Peyrache et al., 2009 & 2010)

% Get function-specific settings
dorandperm      = ft_getopt(cfg, 'dorandperm', false);
nperms          = ft_getopt(cfg, 'nperms', 1000);

%% Get n x t matrix of spikes for training and testing
% Get the train raster
binsize_train = 0.2; % sec
stride_train = 0.04; % sec
toi_train = 'targon'; % we will derive PCs from times when targets are on the screen
[nspk_z_train, t_train, spk_avg_train, spk_sd_train] = ...
  get_raster(kinematics, spike, binsize_train, stride_train, toi_train);

% Get the test raster
binsize_test = 0.02; % sec
stride_test = 0.004; % sec
toi_test = 'all'; % we will apply the projectors to all time bins
[nspk_z_test, t_test, spk_avg_test, spk_sd_test] = ...
  get_raster(kinematics, spike, binsize_test, stride_test, toi_test);

%% Get the reactivations
react = pca_react(nspk_z_train, nspk_z_test);
react.t_test = t_test;
react.input.train.nspk_avg = spk_avg_train;
react.input.train.nspk_sd = spk_sd_train;
react.input.test.nspk_avg = spk_avg_test;
react.input.test.nspk_sd = spk_sd_test;

%% Get the reactivation rate for each PC, each break
% use the 99th percentile of reactivation strengths as the threshold 
% for detecting reactivation events
p_thresh = 0.99; 

[rr_pc, r_thresh] = get_reactrate(react, kinematics, p_thresh);

react.rate_pc_brk = rr_pc;

%% DO THE RANDOM PERMUTATIONS
if dorandperm
  rr_pc_rp = nan(nperms, size(rr_pc, 1), size(rr_pc, 2));
  
  % first get all of the spike times
  ts_all = [spike.ts];
  
  for n = 1:nperms
    % randomly assign new unit identities to each spike
    unit_id_rp = ceil(length(spike)*rand(1, length(ts_all)));
    
    % convert the new spike identities back into a spike structure
    spike_rp = [];
    for i = 1:length(spike)
      spike_rp(i).ts = ts_all(find(unit_id_rp == i));
    end
    
    [nspk_z_rp, t_rp] = ...
      get_raster(kinematics, spike_rp, binsize_test, stride_test, toi_test);
    
    react_rp = pca_react(nspk_z_train, nspk_z_rp);
    react_rp.t_test = t_rp;
    
    rr_pc_rp(n, :, :) = get_reactrate(react_rp, kinematics, [], r_thresh);
    
    react_rp.cfg = cfg;
  end
  
  react.rate_pc_brk_rp = rr_pc_rp;
end

%% SUBFUNCTIONS
function [nspk_z, t_nspk, spk_avg, spk_sd] = get_raster(kinematics, spike, binsize, stride, toi)

% Create a spike count matrix of the entire recording period
[n_spk_all, t_all] = make_raster_stride(spike, binsize, stride);

% zscore the spike count matrix
spk_avg = nanmean(n_spk_all, 2);
spk_sd = nanstd(n_spk_all, 0, 2);

n_spk_all_z = zscore(n_spk_all, 0, 2);

% Get the relevant timepoints
cfg_toi = [];
cfg_toi.kinematics = kinematics;
is_toi = get_is_epoch(cfg_toi, t_all, toi);
nspk_z = n_spk_all_z(:, is_toi);

t_nspk = t_all(is_toi);


function [rr_pc, r_thresh] = get_reactrate(react, kinematics, p_thresh, r_thresh)

pc = find(react.pca.sig_pcs);

if nargin == 4
  get_new_thresh = false;
else
  get_new_thresh = true;
  r_thresh = nan(1, length(pc));
end

rr_pc = nan(length(pc), length(kinematics.break));

cfg_toi = [];
cfg_toi.kinematics = kinematics;

is_brk = get_is_epoch(cfg_toi, react.t_test, 'break');
is_brkmv = get_is_epoch(cfg_toi, react.t_test, 'brkmove');

for p = 1:length(pc)
  if get_new_thresh
    react_all = react.react(pc(p), is_brk & ~is_brkmv);
    r_sort = sort(react_all, 'ascend');
    r_thresh(p) = r_sort(round(p_thresh*length(r_sort)));
  end
  
  for k = 1:length(kinematics.break)
    dur_mv = sum(kinematics.break(k).t_move_end - kinematics.break(k).t_move_start);
    dur_tot = kinematics.break(k).t_end - kinematics.break(k).t_start;
    dur_nomv = dur_tot - dur_mv;
    
    is_brk_k = get_is_epoch(cfg_toi, react.t_test, ['break' num2str(k)]);
    
    n_react = sum(react.react(pc(p), :) > r_thresh(p) & is_brk_k & ~is_brkmv);
    rr_pc(p, k) = n_react/dur_nomv;
  end
end






