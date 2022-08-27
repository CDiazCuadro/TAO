cvx_setup

load A.asc
load b.asc
lambda=1;

%planteo del problema dentro del entorno cvx
cvx_begin 
cvx_precision high
variables x(2) t(2)  
minimize(0.5*(A*x-b)'*(A*x-b)+lambda*sum(t))
subject to
    t>=x 
    t>=-x
cvx_end
opt_val=0.5*(A*x-b)'*(A*x-b)+lambda*sum(t)
x
t
