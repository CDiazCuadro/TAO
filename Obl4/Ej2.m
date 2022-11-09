%% EJERCICIO 2 - OBLIGATORIO 4
%
% Proximal Gradient Method
%
% x(next) = prox[lambda,g](x(actual) - lambda.Grad_f(x(actual)) )
%%
close all, clear all
% Loading data
load A.asc
load b.asc

% Initial parameters
lambda  = 0.15      ;
tol_dif = 1e-4      ;
dif     = tol_dif +1; 
max_it  = 1000000   ;
alpha = 1/norm(A'*A);
al    = alpha*lambda;

% What you wanna run?
RUN_PGD  = 1;
RUN_ADMM = 0;

% General Initializing 
n = size(A,2)     ;
rand ("seed", 10) ; % Para que x0 sea siempre la misma
x0 = rand(n,1)    ;
F0 = 0.5*norm(A*x0 - b)^2 + lambda*norm(x0,1);
x_next = zeros(2,1);

if RUN_PGD
 % PGD Initializing
  X_PGD  = [x0] ;
  it_PGD = 0    ;
  F_PGD = [F0]  ;
  DIF_PGD = []  ;
  tic()
  while tol_dif < dif && it_PGD < max_it
    
    % Defining actual variables
    x_act      = X_PGD(:,end)                 ;
    grad_f_act = A'*(A*x_act - b)             ;
    F_act      = 0.5*norm(A*x_act - b)^2 +...
                 lambda*norm(x_act,1)         ; 

    % Computing
    x_aux  = x_act - alpha*grad_f_act ;
    for i = 1:n
      if x_aux(i) > al
        x_next(i) = x_aux(i) - al;
      elseif  x_aux(i) < -al
        x_next(i) = x_aux(i) + al;
      else
        x_next(i) = 0;
      endif
    endfor
    F_next = 0.5*norm(A*x_next - b)^2 + lambda*norm(x_next,1);
    dif    = abs(F_next - F_act) ;
    DIF_PGD = [DIF_PGD,dif];

    % Storaging & Updating
    X_PGD= [X_PGD, x_next];
    F_PGD = [F_PGD, F_next];
    it_PGD += 1        ;
  endwhile
  time_PDG = toc();

  % Ploting 
  figure
  plot(F_PGD,'o-', "markersize", 3)
  xlabel('Numero de iteracion')
  ylabel('Función objetivo para x(it)')
  grid minor on
  
endif

if RUN_ADMM
 % PGD Initializing
  X_ADMM  = [x0] ;
  it_ADMM = 0    ;
  F_ADMM = [F0]  ;
  DIF_ADMM = []  ;
  tic()
  while tol_dif < dif && it_ADMM < max_it
    
    % Defining actual variables

    % Computing

    % Storaging & Updating
  
  endwhile
  time_ADMM = toc();
  % Ploting 
  
endif