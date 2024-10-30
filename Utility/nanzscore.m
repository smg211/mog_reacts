function dat_out = nanzscore(dat_in, dim)
% zscores data ignoring nan values
% use as: dat_out = nanzscore(dat_in, dim)

dat_out = (dat_in - nanmean(dat_in, dim))./nanstd(dat_in, 0, dim);