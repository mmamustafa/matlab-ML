function t=evaluate_poly(x, coeff)
  sys_ord = length(coeff) - 1;
  t = coeff(end) * ones(size(x));
  for i=sys_ord:-1:1
    t = t + coeff(sys_ord-i+1) * x.^i;
  end
  