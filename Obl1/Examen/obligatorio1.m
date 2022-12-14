%% Ejercicio 4:
clc
close all
clear all

A= [-4.100000000000000000e+01 2.000000000000000000e+01;...
-4.600000000000000000e+01 -8.000000000000000000e+00;...
-5.000000000000000000e+00 -3.300000000000000000e+01;...
-5.500000000000000000e+01 1.000000000000000000e+00;...
-5.500000000000000000e+01 -6.000000000000000000e+00];

b= [8.000000000000000000e+00
5.000000000000000000e+00
-3.000000000000000000e+00
1.000000000000000000e+01
4.000000000000000000e+00];

x_min = A\b;

%% Parte c - Descenso por gradiente:
%
% Datos genericos para todas las partes:
tol_f = 1e-5                ; % Tolerancia para f(x)
it_max = 5000               ; % Determino el maximo de iteraciones
x_0 = 100*rand(2,1)         ; % Punto de partida del algoritmo
grad_0 = (2*(A*x_0 - b)'*A)'; % Calculo el gradiente en x_0
f_0 = (norm(A*x_0 -b))^2    ; % calculo el valor funcional de x_0
var_0 = tol_f + 1           ; % Aseguro que var > tol la primera vez
rel_error_0 = norm(x_min - x_0)/norm(x_min) ; 

% --> a) paso fijo   s=1/||A||?_2
%
% Variables que dependen de cada caso:
s_a = 1/(2*norm(A)^2)         ; % Defino el paso
it_a = 0                      ; % Inicializo las iteraciones
x_a_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
grad_a_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
var_a_hist = []               ; % inicializo la matriz de historial de variacion
var_a = var_0                 ; % Aseguro que var > tol la primera vez
rel_error_a = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
iterations_a = [0]            ;

while var_a > tol_f && it_a < it_max
  it_a = it_a + 1                                 ; % Actualizo el numero de iteracion
  d_a       = -grad_a_hist(:,end)           ; % Defino direccion opuesto al gradiente
  x_a_act   = x_a_hist(:,end)               ; % Obtengo la posicion actual   
  x_a_next  = x_a_act + s_a*d_a             ; % Calculo la siguiente posicion
  grad_a_next = (2*(A*x_a_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
  % Almaceno valores  
  x_a_hist = [x_a_hist, x_a_next]           ;
  grad_a_hist = [grad_a_hist , grad_a_next] ;
  % Calculo variacion
  f_a_act = (norm(A*x_a_act -b))^2          ;
  f_a_next = (norm(A*x_a_next -b))^2        ;
  var_a = norm(f_a_next-f_a_act)/norm(f_0)  ;
  var_a_hist = [var_a_hist, var_a]          ;
  error_a_it = norm(x_min-x_a_next)/norm(x_min) ;
  rel_error_a = [rel_error_a,error_a_it]        ;
  iterations_a = [iterations_a,it_a]            ;

end

%
% --> b) paso decreciente   s=0.001*1/k
%
% Variables que dependen de cada caso:
it_b = 0                      ; % Inicializo las iteraciones
x_b_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
grad_b_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
var_b_hist = []               ; % inicializo la matriz de historial de variacion
var_b = var_0                 ; % Aseguro que var > tol la primera vez
rel_error_b = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
iterations_b = [0]            ;


while var_b > tol_f && it_b < it_max
  it_b = it_b + 1                                 ; % Actualizo el numero de iteracion
  s_b = 0.001/it_b                          ; % Calculo el paso
  d_b       = -grad_b_hist(:,end)           ; % Defino direccion opuesto al gradiente
  x_b_act   = x_b_hist(:,end)               ; % Obtengo la posicion actual   
  x_b_next  = x_b_act + s_b*d_b             ; % Calculo la siguiente posicion
  grad_b_next = (2*(A*x_b_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
  % Almaceno valores  
  x_b_hist = [x_b_hist, x_b_next]           ;
  grad_b_hist = [grad_b_hist , grad_b_next] ;
  % Calculo variacion
  f_b_act = (norm(A*x_b_act -b))^2          ;
  f_b_next = (norm(A*x_b_next -b))^2        ;
  var_b = norm(f_b_next-f_b_act)/norm(f_0)  ;
  var_b_hist = [var_b_hist, var_b]          ;
  error_b_it = norm(x_min-x_b_next)/norm(x_min) ;
  rel_error_b = [rel_error_b,error_b_it]        ;
  iterations_b = [iterations_b,it_b]            ;
end


%
% --> c) paso exacto (line search)   
%      - Se que x*=A\b  --> busco s: xk +sk*dk = x*
%                       --> sk= (A\b-xk)/dk
%
%
% Variables que dependen de cada caso:
it_c = 0                      ; % Inicializo las iteraciones
x_c_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
grad_c_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
var_c_hist = []               ; % inicializo la matriz de historial de variacion
var_c = var_0                 ; % Aseguro que var > tol la primera vez
rel_error_c = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
iterations_c = [0]            ;

while var_c > tol_f && it_c < it_max
  it_c = it_c + 1                                 ; % Actualizo el numero de iteracion
  d_c       = -grad_c_hist(:,end)           ; % Defino direccion opuesto al gradiente
  x_c_act   = x_c_hist(:,end)               ; % Obtengo la posicion actual   
  s_c = (A\b - x_c_act)/d_c                 ; % Calculo el paso
  x_c_next  = x_c_act + s_c*d_c             ; % Calculo la siguiente posicion
  grad_c_next = (2*(A*x_c_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
  % Almaceno valores  
  x_c_hist = [x_c_hist, x_c_next]           ;
  grad_c_hist = [grad_c_hist , grad_c_next] ;
  % Calculo variacion
  f_c_act = (norm(A*x_c_act -b))^2          ;
  f_c_next = (norm(A*x_c_next -b))^2        ;
  var_c = norm(f_c_next-f_c_act)/norm(f_0)  ;
  var_c_hist = [var_c_hist, var_c]          ;
  error_c_it = norm(x_min-x_c_next)/norm(x_min) ;
  rel_error_c = [rel_error_c,error_c_it]        ;
  iterations_c = [iterations_c,it_c]            ;
end


%
% --> d) regla de Armijo    
%      - Defino: sig = 0.1
%                beta = 1/5       
%      --> iterar j hasta que: f_old - f_next > - sig*b^j*s*gradF*dk
%
% Variables que dependen de cada caso:
it_d = 0                      ; % Inicializo las iteraciones
x_d_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
grad_d_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
var_d_hist = []               ; % inicializo la matriz de historial de variacion
var_d = var_0                 ; % Aseguro que var > tol la primera vez
sig = 0.1                     ; % Parametro 1 de armijo
beta = 1/5                    ; % Parametro 2 de armijo
s_d = 1                       ; % Paso original
rel_error_d = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
iterations_d = [0]            ;


while var_d > tol_f && it_d < it_max
  it_d =it_d + 1                                 ; % Actualizo el numero de iteracion
  d_d       = -grad_d_hist(:,end)           ; % Defino direccion opuesto al gradiente
  x_d_act   = x_d_hist(:,end)               ; % Obtengo la posicion actual 
  j = 1                                     ; % iterador de armijo
  gan = -sig*beta^(j)*s_d*d_d'*d_d          ; % Ganancia si fuera lineal
  x_d_next_j = x_d_act + beta^(j)*s_d*d_d   ; % Siguiente posicion en la iteracion j=1
  f_d_act = (norm(A*x_d_act -b))^2          ; % Valor funcional de la posicion actual
  f_d_next_j = (norm(A*x_d_next_j -b))^2    ; % Valor funcional del primer intento de armijo
  gan_real = f_d_act - f_d_next_j           ;
  while  gan_real < gan
    j = j+ 1                                  ; % Aumento el iterador de armijo
    gan = -sig*beta^(j)*s_d*d_d'*d_d        ; % Re-calculo la ganancia
    x_d_next_j = x_d_act + beta^(j)*s_d*d_d ; % Siguiente posicion en la iteracion
    f_d_next_j = (norm(A*x_d_next_j -b))^2  ; % Valor funcional del primer intento de armijo
    gan_real = f_d_act - f_d_next_j         ;
  end 
  x_d_next  = x_d_next_j                    ; % Defino la posicion
  grad_d_next = (2*(A*x_d_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
  % Almaceno valores  
  x_d_hist = [x_d_hist, x_d_next]           ;
  grad_d_hist = [grad_d_hist , grad_d_next] ;
  % Calculo variacion
  f_d_next = f_d_next_j                     ;
  var_d = norm(f_d_next-f_d_act)/norm(f_0)  ;
  var_d_hist = [var_d_hist, var_d]          ;
  error_d_it = norm(x_min-x_d_next)/norm(x_min) ;
  rel_error_d = [rel_error_d,error_d_it]        ;
  iterations_d = [iterations_d,it_d]            ;
  end

figure
hold on
plot(iterations_a,rel_error_a,'*-')
plot(iterations_c,rel_error_c,'*-')
plot(iterations_d,rel_error_d,'*-')
xl = xlabel('Numero de iteracion')
yl = ylabel('Error relativo')
lg = legend('Paso fijo','Paso exacto','Regla de armijo')
grid minor
set (xl, "fontsize", 16);
set (yl, "fontsize", 16);
set (lg, "fontsize", 16);

figure
plot(iterations_b,rel_error_b,'*-')
xl = xlabel('Numero de iteracion')
yl = ylabel('Error relativo')
lg = legend('Paso decreciente')
grid minor
set (xl, "fontsize", 16);
set (yl, "fontsize", 16);
set (lg, "fontsize", 16);
