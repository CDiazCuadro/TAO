%%%%% Funcion para calcular el gradiente de h(ê)
function grad = grad_h(e)
   n = length(e);
   grad = zeros(n,1);
   for idx = 1 : n
       ei = e(idx);
       if ei >= 0
           grad(idx) = (ei^2 + 2*ei)/(ei+1)^2;
       else
           grad(idx) = (3*ei^2+2*ei)/(1-ei)^2;
       end                
   end
end