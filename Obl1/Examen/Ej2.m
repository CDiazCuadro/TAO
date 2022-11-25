%% EJERCICIO 2 - EXAMEN
%
% Proximal Gradient Method
%
% x(next) = prox[norm_inf](x(actual) - alpha.Grad_f(x(actual)) )
%%
close all, clear all, clc
% Loading data
% load A_May.txt 
% A = A_May;
% load b.txt
% % b=b*0;
load x_antilasso.txt
b=[2;3];
A=[5 1; -1 0];
% x_antilasso = inv(A)*b;

 


%
% Preparing figures
figure (1)
hold on
grid  minor
xlabel('Numero de iteracion')
ylabel('Funcion')
%
figure (2)
hold on
grid  minor
xlabel('Numero de iteracion')
ylabel('Diferencia de la función de costo')
%
% figure (3)
% hold on
% grid  minor
% xlabel('Numero de iteracion')
% ylabel('Diferencia con la solución X^*')
%
valores=[10];
%
for div = valores
    % Initial parameters
    tol_dif = 1e-4      ;
    dif     = tol_dif +1; 
    max_it  = 10000   ;
    alpha = 0.01         ;
    beta  = 1           ;
    lam   = 1 ;
    %
    % General Initializing 
    n = size(A,2)     ;
    rand ("seed", pi()) ; % Para que x0 sea siempre la misma
    x0 = 100*rand(n,1)    ;
    F0 = (1/2/beta)*norm(A*x0 - b)^2 + norm(x0,Inf);
    % 
    % PGD Initializing
    X  = [x0] ;
    it = 0    ;
    F = [F0]  ;
    DIF = [F0]  ;
%     DIF_SOL=[norm(x_antilasso-x0)];
    DIF_SOL=[norm(x0)];

    %
    while tol_dif < dif && it < max_it
        %
        % Defining actual variables
        x_act      = X(:,end);
        F_act      = (1/2/beta)*norm(A*x_act - b)^2 + norm(x_act,Inf); 
        grad_f_act = (1/beta)*A'*(A*x_act - b);
        % 
        % Computing
        x_aux  = x_act - alpha*grad_f_act;
        x_next = prox_inf(x_aux, lam);
        F_next = (1/2/beta)*norm(A*x_next - b)^2 + norm(x_next,Inf);
        dif    = abs(F_next - F_act) ;
%         dif_sol= norm(x_antilasso - x_next);
        dif_sol= norm(x_next);
        DIF    = [DIF, dif];
        DIF_SOL= [DIF_SOL, dif_sol];
        % 
        % Storaging & Updating
        X = [X, x_next];
        F = [F, F_next];
        it= it + 1     ;
    end
    if it == max_it
        text = ['Con alpha = 1/AtA/',num2str(div),' cortó por max_iter...\n'];
        fprintf(text)
    else
        text = ['Con alpha = 1/AtA/',num2str(div),' cortó por convergencia con ',num2str(it),' iteraciones \n'];
        fprintf(text)
    end
    %
    absisas = [0:length(F)-1];
    figure (1)
    plot(absisas, F)
    figure (2)
    semilogy(absisas, DIF)
%     figure (3)
%     plot(absisas, DIF_SOL)
end
