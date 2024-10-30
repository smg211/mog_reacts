function [p_x, x_pmax] = do_bayes_decode(fxmatrix, n_spk, binsize_sec)

n_t = size(n_spk, 2);
n_x = size(fxmatrix, 2); % number of conditions to decode between
n_unit = size(fxmatrix, 1);

p_x = nan(n_x, n_t);
x_pmax = nan(1, n_t);

for t = 1:n_t
  endprob = nan(1, n_x);
  
  for x = 1:n_x
    productme = 0;
    expme = 0;
    
    for n = 1:n_unit
      ni = n_spk(n, t);
      fx = (fxmatrix(n, x));  % the rate for cell u during condition X
      
      
      if fx ~= 0
        productme = productme + ni*log(fx);
      else
        fx = eps;
        productme = productme + ni*log(fx);
      end
      
      productme = (productme + ni*log(fx));
      expme = (expme) + (fx);
    end % end unit loop
    
    endprob(x) = (productme) + (-binsize_sec.*expme);
  end
  
  [~, i_max] = (max(endprob));
  
  nums = isfinite(endprob);
  nums = find(nums == 1);
  endprob = endprob(nums);
  
  % why is this 12?
  mp = max(endprob(:))-12;
  
  endprob = exp(endprob-mp);
  conv = 1./sum(endprob(~isnan(endprob)), 'all');
  endprob = endprob*conv;
  
  p_x(:, t) = endprob;
  
  x_pmax(t) = i_max;
end














