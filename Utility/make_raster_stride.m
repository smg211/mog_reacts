function [raster, t_raster, t_binedges] = make_raster_stride(spike, binsize, stride)
% Use as:
%     [raster, t_raster] = make_raster_stride(spike, binsize, stride, t_start, t_end)
%
% Where:
%     spike   = structure, with N_unit elements where each element contains a
%             field named ts which contains all spike times for that unit (in sec)
%     binsize = scalar, width of bins (seconds)
%     stride  = scalar, stride length for binning (seconds)
% 
% Sandon Griffin (2024)

t_start = 0;
t_end = max([spike.ts]);

% get the binedges using the stride length
binedges = t_start:stride:t_end+stride;

% rebin the events
dat_out = [];
t_binedges = [];
for u = 1:length(spike)
  hc_stride = histcounts(spike(u).ts, binedges);
  hc_window = zeros(1, length(binedges)-1);
  for i = 1:(binsize/stride)
    hc_window = hc_window + [hc_stride(i:end), zeros(1, i-1)];
  end
  raster(u, :) = hc_window(1:end-i);
end

t_raster = (t_start+binsize/2):stride:(t_end-binsize/2);
t_binedges(:, 1) = t_raster-binsize/2;
t_binedges(:, 2) = t_raster+binsize/2;