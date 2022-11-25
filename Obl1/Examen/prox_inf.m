function Xprox = prox_inf(x, lambda)
    Xprox = x - lambda*proyL1(x./lambda);
end