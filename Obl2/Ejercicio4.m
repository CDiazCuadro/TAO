%% Ejercicio 4

clear all, close all, clc
% Carga de datos
load A.txt 
load B.txt
load C.txt 

load alpha.txt 
load beta.txt 
load gamma.txt

%% Defino matrices auxiliares D, delta y H
n = length(A);

% Matriz D
D = zeros (3*n,3*n);
D(1:n,1:n)              = A;
D(n+1:2*n,n+1:2*n)      = B;
D(2*n+1:3*n,2*n+1:3*n ) = C;

% Matriz H 
H = zeros(2*n,3*n);
H(1:n,1:n)            = eye(n);
H(1:n,n+1:2*n)        = -eye(n);
H(n+1:2*n, 1:n)       = eye(n);
H(n+1:2*n, 2*n+1:3*n) = -eye(n);

% Vector delta
delta = zeros(3*n,1);
delta = [alpha; beta; gamma];  

%% Calculo solucion exacta x_op = inv(M)*m
% Defino M, m y mu
M  = A'*A + B'*B + C'*C;
m  = A'*alpha + B'*beta + C'*gamma;
mu = alpha'*alpha + beta'*beta + gamma'*gamma;
x_op = inv(M)*m ;

%% Calculo los multiplicadores lambda_op
A = H*(D'*D)^(-1)*D'*delta;
B = H*(D'*D)^(-1)*H';
lambda_op = inv(B)*A;

%% Solucion Nummerica:
diff_max = 1e-8;
k_max = 100 ;
t0 = 1/norm(D);

%%
% Parte E) lambda = lambda_op
%          tk = t02^k.

E_lambda = lambda_op  ;
E_t0 = t0             ;
E_t_hist = [E_t0]     ; 
E_k      = 0          ;  
E_diff  = diff_max + 1;
E_w0 = inv( (1/2)*(D'*D + (D'*D)') + (E_t0/2)*(H'*H + (H'*H)') )*(D'*delta - H'*E_lambda) ;
E_w_hist  = [E_w0]      ;
E_x0 = (E_w0(1:10) + E_w0(11:20) + E_w0(21:30))/3;
E_e0 = norm(E_x0-x_op)/norm(x_op);
E_e_hist  = [E_e0]      ;     
E_diff_hist = [E_diff];

while E_diff > diff_max && E_k < k_max
      E_k = E_k + 1 ;
      E_t = E_t0*2^E_k;
      w_old = E_w_hist(:,end);
      w = inv( (1/2)*(D'*D + (D'*D)') + (E_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*E_lambda) ;
      E_diff = norm(w - w_old)/norm(w);
      E_diff_hist = [E_diff_hist,E_diff]  ;
      E_w_hist = [E_w_hist,w];
      x_num = (w(1:10) + w(11:20) + w(21:30))/3 ;
      e_act = norm(x_num-x_op)/norm(x_op);
      E_e_hist= [E_e_hist, e_act] ;
      E_t_hist = [E_t_hist, E_t];

end
E_x = [0:E_k];

%%
% Parte F) lambda = lambda_op*0
%          tk = t02^k.

F_lambda = lambda_op*0;
F_t0 = t0             ;
F_t_hist = [F_t0]     ; 
F_k      = 0          ;  
F_diff  = diff_max + 1;
F_w0 = inv( (1/2)*(D'*D + (D'*D)') + (F_t0/2)*(H'*H + (H'*H)') )*(D'*delta - H'*F_lambda) ;
F_w_hist  = [F_w0]      ;
F_x0 = (F_w0(1:10) + F_w0(11:20) + F_w0(21:30))/3;
F_e0 = norm(F_x0-x_op)/norm(x_op);
F_e_hist  = [F_e0]      ;     
F_diff_hist = [F_diff]  ;

while F_diff > diff_max && F_k < k_max
      F_k = F_k + 1 ;
      F_t = F_t0*2^F_k;
      w_old = F_w_hist(:,end);
      w = inv( (1/2)*(D'*D + (D'*D)') + (F_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*F_lambda) ;
      F_diff = norm(w - w_old)/norm(w);
      F_diff_hist = [F_diff_hist,F_diff]  ;
      F_w_hist = [F_w_hist,w];
      x_num = (w(1:10) + w(11:20) + w(21:30))/3 ;
      e_act = norm(x_num-x_op)/norm(x_op);
      F_e_hist= [F_e_hist, e_act] ;
      F_t_hist = [F_t_hist, F_t];
end 
F_x = [0:F_k];

%%
% Parte G1) lambda0 = lambda_op*0
%          lambdak = lambda(k-1)+tHw
%          t = 10*t0

G1_lambda0 = lambda_op*0     ;
G1_lambda_hist = [G1_lambda0];
G1_t0 = t0                   ;
G1_k = 0                     ;
G1_t = 10*G1_t0              ;
G1_diff  = diff_max + 1      ;
G1_w0 = inv( (1/2)*(D'*D + (D'*D)') + (G1_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*G1_lambda0);
G1_w_hist  = [G1_w0]            ;
G1_x0 = (G1_w0(1:10) + G1_w0(11:20) + G1_w0(21:30))/3;
G1_e0 = norm(G1_x0-x_op)/norm(x_op);
G1_e_hist  = [G1_e0]            ;     
G1_diff_hist = [G1_diff]     ;

while G1_diff > diff_max && G1_k < k_max
      G1_k = G1_k + 1 ;
      w_old = G1_w_hist(:,end);
      G1_lambda = G1_lambda_hist(:,end) + G1_t*H*w_old;
      w = inv( (1/2)*(D'*D + (D'*D)') + (G1_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*G1_lambda);
      G1_diff = norm(w - w_old)/norm(w);
      x_num = (w(1:10) + w(11:20) + w(21:30))/3; 
      e_act = norm(x_num-x_op)/norm(x_op);
      
      % Almaceno
      G1_lambda_hist = [G1_lambda_hist, G1_lambda] ; 
      G1_diff_hist = [G1_diff_hist,G1_diff]        ;
      G1_w_hist = [G1_w_hist,w]                   ;
      G1_e_hist= [G1_e_hist, e_act]               ;
end 
G1_x = [0:G1_k];

%%
% Parte G2) lambda0 = lambda_op*0
%          lambdak = lambda(k-1)+tHw
%          t = 1000*t0

G2_lambda0 = lambda_op*0     ;
G2_lambda_hist = [G2_lambda0];
G2_t0 = t0                   ;
G2_k = 0                     ;
G2_t = 1000*G2_t0            ;
G2_diff  = diff_max + 1      ;
G2_w0 = inv( (1/2)*(D'*D +(D'*D)') + (G2_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*G2_lambda0);
G2_w_hist  = [G2_w0]            ;
G2_x0 = (G2_w0(1:10) + G2_w0(11:20) + G2_w0(21:30))/3;
G2_e0 = norm(G2_x0-x_op)/norm(x_op);
G2_e_hist  = [G2_e0]            ;     
G2_diff_hist = [G2_diff]     ;

while G2_diff > diff_max && G2_k < k_max
      G2_k = G2_k + 1 ;
      w_old = G2_w_hist(:,end);
      G2_lambda = G2_lambda_hist(:,end) + G2_t*H*w_old;
      w = inv( (1/2)*(D'*D +(D'*D)') + (G2_t/2)*(H'*H + (H'*H)'))*(D'*delta - H'*G2_lambda);
      G2_diff = norm(w - w_old)/norm(w);
      x_num = (w(1:10) + w(11:20) + w(21:30))/3; 
      e_act = norm(x_num-x_op)/norm(x_op);
      
      % Almaceno
      G2_lambda_hist = [G2_lambda_hist, G2_lambda] ; 
      G2_diff_hist = [G2_diff_hist,G2_diff]        ;
      G2_w_hist = [G2_w_hist,w]                    ;
      G2_e_hist= [G2_e_hist, e_act]                ;
end 
G2_x = [0:G2_k];

%%
% Parte H) lambda0 = lambda_op*0
%          lambdak = lambda(k-1)+tHw
%          t0 = 1000*t0
%          t = t02^k.

H_lambda0 = lambda_op*0     ;
H_lambda_hist = [H_lambda0];
H_t0 = t0;
H_t_hist = [H_t0];
H_k = 0                     ;
H_diff  = diff_max + 1      ;
H_w0 = inv( (1/2)*(D'*D + (D'*D)') + (H_t0/2)*(H'*H + (H'*H)') )*(D'*delta - H'*H_lambda0);
H_w_hist  = [H_w0]            ;
H_x0 = (H_w0(1:10) + H_w0(11:20) + H_w0(21:30))/3;
H_e0 = norm(H_x0-x_op)/norm(x_op);
H_e_hist  = [H_e0]            ;     
H_diff_hist = [H_diff]     ;

while H_diff > diff_max && H_k < k_max
      H_k = H_k + 1 ;
      w_old = H_w_hist(:,end);
      H_t = H_t0*2^H_k;
      H_lambda = H_lambda_hist(:,end) + H_t*H*w_old;
      w = inv( (1/2)*(D'*D + (D'*D)') + (H_t/2)*(H'*H + (H'*H)') )*(D'*delta - H'*H_lambda);
      H_diff = norm(w - w_old)/norm(w);
      x_num = (w(1:10) + w(11:20) + w(21:30))/3; 
      e_act = norm(x_num-x_op)/norm(x_op);
      
      % Almaceno
      H_lambda_hist = [H_lambda_hist, H_lambda] ; 
      H_diff_hist = [H_diff_hist,H_diff]        ;
      H_w_hist = [H_w_hist,w]                   ;
      H_e_hist= [H_e_hist, e_act]               ;
end 
H_x = [0:H_k];

%%

figure (1)
semilogy(E_x , E_e_hist ,  'o-',...
         F_x , F_e_hist ,  'o-',...
         G1_x, G1_e_hist,  'o-',...
	     G2_x, G2_e_hist,  'o-',...
         H_x , H_e_hist ,  'o-')

legend('Parte e) {\lambda} = {\lambda}* || {\tau}^k={\tau}_02^k ',...
            'Parte f) {\lambda} = 0 || {\tau}^k={\tau}_02^k ',...
            'Parte g1) {\lambda}^k = {\lambda}^{k-1}+{\tau}H{\omega} || {\tau}^k=10{\tau}_0 ',...
            'Parte g2) {\lambda}^k = {\lambda}^{k-1}+{\tau}H{\omega} || {\tau}^k=1000{\tau}_0 ',...
            'Parte h) {\lambda}^k = {\lambda}^{k-1}+{\tau}H{\omega} || {\tau}^k={\tau}_02^k ')

xlabel('Numero de iteracion')
ylabel('Error relativo ||x^k-x*||/||x*||')
grid minor