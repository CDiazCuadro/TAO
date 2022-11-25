function[u] = shrink(v, lambdaIn)


if strcmp(class(v), 'gpuArray')
  lambda = gpuArray(lambdaIn);
else
  lambda = lambdaIn;
end

if isscalar(lambda),
  u = sign(v).*max(0, abs(v) - lambda);
else
  u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
end

  
return
