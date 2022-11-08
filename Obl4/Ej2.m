%% EJERCICIO 2 - OBLIGATORIO 4
%
% Proximal Gradient Method
%
% x(next) = prox[lamda,g](x(actual) - lambda.Grad_f(x(actual)) )
%%
% Loading data
load A.asc
load b.asc

% Initial parameters
lambda  = 0.15;
tol_dif = 1e-4;
max_it  = 100 ;

% Initializing 
rng('default') % Para que x0 sea siempre la misma
x0 = rand(size(A,2),1)  ;
X  = [x0]               ;
it = 0                  ;


alpha = 1/norm(A'*A);
grad_f = 
x_aux = X(end) - alpha*grad_f