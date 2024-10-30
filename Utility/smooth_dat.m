function dat_sm = smooth_dat(dat, sm_kernel, method)
dat_sm = dat;
if strcmp(method, 'moving')
  % deal with padding
  pad = ceil(sm_kernel/2);
  dat_sm = ft_preproc_padding(dat_sm, 'localmean', pad);

  for d = 1:size(dat_sm, 1)
    dat_sm(d, :) = smooth(dat_sm(d, :), sm_kernel, 'moving');
  end

  % cut the eges
  dat_sm = ft_preproc_padding(dat_sm, 'remove', pad);
elseif strcmp(method, 'gaussian')
  dat_sm = smoother(dat, sm_kernel, 1);
end
