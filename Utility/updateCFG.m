function cfg_out = updateCFG(cfg_in, vars)

cfg_out = cfg_in;
for i = 1:length(vars)
  cfg_out.(vars{i}) = evalin('caller', vars{i});
end