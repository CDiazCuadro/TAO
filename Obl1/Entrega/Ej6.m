%% Ejercicio 6

clc
close all
clear all

f = @(x,y) 5*x*x + 5*y*y +5*x - 3*y - 6*x*y +5/4;
df_x = @(x,y) 10*x + 5 -6*y;
df_y = @(x,y) 10*y -3 -6*x ;

% Datos genericos para todas las partes:
x_min = [-1/2 ; 0]                          ; % Solucion
R = 1/4                                     ; % Radio de la region factible 
tol_var = 1e-5                                ; % Tolerancia para (x_next-x_act)
tol_real = 1e-5                             ; % Tolerancia con el error en la solucion real
it_max = 50000                               ; % Determino el maximo de iteraciones
x_0 = R*rand(2,1)                           ; % Punto de partida del algoritmo
grad_0 = [df_x(x_0(1),x_0(2)) ;...
          df_y(x_0(1),x_0(2))]              ; % Calculo el gradiente en x_0
f_0 = f(x_0(1),x_0(2))                      ; % calculo el valor funcional de x_0
var_0 = tol_var + 1                           ; % Aseguro que var > tol la primera vez
rel_error_0 = norm(x_min - x_0)/norm(x_min) ; 


%
% --> b) paso decreciente   s=1/k
%                           alpha = 1
% Variables que dependen de cada caso:
it_b = 0                        ; % Inicializo las iteraciones
x_b_hist = [x_0]                ; % inicializo la matriz de historial de posiciones
grad_b_hist = [grad_0]          ; % inicializo la matriz de almacenamiento de gradientes
xt_b_0 = proyect(x_0-grad_0 ,R) ; % Posicion inicial hacia donde apunto, proyectada.
xt_b_hist = [xt_b_0]            ; % % inicializo el vector de historial de posiciones hacia donde apunto
var_b_hist = []                 ; % inicializo la matriz de historial de variacion
var_b = var_0                   ; % Aseguro que var > tol la primera vez
rel_error_b = [rel_error_0]     ; % inicializo el vector de historial de errores relativos
iterations_b = [0]              ;
error_b_it = rel_error_0        ;
f_b_hist = [f_0]                ; 

while error_b_it > tol_real && it_b < it_max && var_b > tol_var 
  it_b += 1                                     ; % Actualizo el numero de iteracion
  s_b = 1/it_b                                  ; % Calculo el paso
  x_b_act   = x_b_hist(:,end)                   ; % Obtengo la posicion actual   
  grad_b_act = grad_b_hist(:,end)               ; 
  x_b_try = x_b_act - s_b*grad_b_act            ;
  x_b_next  = proyect(x_b_try,R);          ; % Calculo la siguiente posicion
  grad_b_next = [df_x(x_b_next(1),x_b_next(2))  ;...
                 df_y(x_b_next(1),x_b_next(2))] ; % Calculo el nuevo cetor de gradiente
  % Almaceno valores  
  x_b_hist = [x_b_hist, x_b_next]               ;
  grad_b_hist = [grad_b_hist , grad_b_next]     ;
  % Calculo variacion
  var_b = norm(x_b_next-x_b_act)  ;
  var_b_hist = [var_b_hist, var_b]          ;
  % Calculo errores
  error_b_it = norm(x_min-x_b_next) ;
  rel_error_b = [rel_error_b,error_b_it]        ;
  % Cuento it
  iterations_b = [iterations_b,it_b]            ;
  % Calculo funcional
  funcional_b = f(x_b_next(1),x_b_next(2))      ;
  f_b_hist = [f_b_hist, funcional_b]  ;
end


%%
%% --> c) paso exacto (line search)   
%%      - Se que x*=A\b  --> busco s: xk +sk*dk = x*
%%                       --> sk= (A\b-xk)/dk
%%
%%
%% Variables que dependen de cada caso:
%it_c = 0                      ; % Inicializo las iteraciones
%x_c_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
%grad_c_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
%var_c_hist = []               ; % inicializo la matriz de historial de variacion
%var_c = var_0                 ; % Aseguro que var > tol la primera vez
%rel_error_c = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
%iterations_c = [0]            ;
%
%while var_c > tol_f && it_c < it_max
%  it_c += 1                                 ; % Actualizo el numero de iteracion
%  d_c       = -grad_c_hist(:,end)           ; % Defino direccion opuesto al gradiente
%  x_c_act   = x_c_hist(:,end)               ; % Obtengo la posicion actual   
%  s_c = (A\b - x_c_act)/d_c                 ; % Calculo el paso
%  x_c_next  = x_c_act + s_c*d_c             ; % Calculo la siguiente posicion
%  grad_c_next = (2*(A*x_c_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
%  % Almaceno valores  
%  x_c_hist = [x_c_hist, x_c_next]           ;
%  grad_c_hist = [grad_c_hist , grad_c_next] ;
%  % Calculo variacion
%  f_c_act = (norm(A*x_c_act -b))^2          ;
%  f_c_next = (norm(A*x_c_next -b))^2        ;
%  var_c = norm(f_c_next-f_c_act)/norm(f_0)  ;
%  var_c_hist = [var_c_hist, var_c]          ;
%  error_c_it = norm(x_min-x_c_next)/norm(x_min) ;
%  rel_error_c = [rel_error_c,error_c_it]        ;
%  iterations_c = [iterations_c,it_c]            ;
%end
%
%
%%
%% --> d) regla de Armijo    
%%      - Defino: sig = 0.1
%%                beta = 1/5       
%%      --> iterar j hasta que: f_old - f_next > - sig*b^j*s*gradF*dk
%%
%% Variables que dependen de cada caso:
%it_d = 0                      ; % Inicializo las iteraciones
%x_d_hist = [x_0]              ; % inicializo la matriz de historial de posiciones
%grad_d_hist = [grad_0]        ; % inicializo la matriz de almacenamiento de gradientes
%var_d_hist = []               ; % inicializo la matriz de historial de variacion
%var_d = var_0                 ; % Aseguro que var > tol la primera vez
%sig = 0.1                     ; % Parametro 1 de armijo
%beta = 1/5                    ; % Parametro 2 de armijo
%s_d = 1                       ; % Paso original
%rel_error_d = [rel_error_0]   ; % inicializo el vector de historial de errores relativos
%iterations_d = [0]            ;
%
%
%while var_d > tol_f && it_d < it_max
%  it_d += 1                                 ; % Actualizo el numero de iteracion
%  d_d       = -grad_d_hist(:,end)           ; % Defino direccion opuesto al gradiente
%  x_d_act   = x_d_hist(:,end)               ; % Obtengo la posicion actual 
%  j = 1                                     ; % iterador de armijo
%  gan = -sig*beta^(j)*s_d*d_d'*d_d          ; % Ganancia si fuera lineal
%  x_d_next_j = x_d_act + beta^(j)*s_d*d_d   ; % Siguiente posicion en la iteracion j=1
%  f_d_act = (norm(A*x_d_act -b))^2          ; % Valor funcional de la posicion actual
%  f_d_next_j = (norm(A*x_d_next_j -b))^2    ; % Valor funcional del primer intento de armijo
%  gan_real = f_d_act - f_d_next_j           ;
%  while  gan_real < gan
%    j = j+ 1                                  ; % Aumento el iterador de armijo
%    gan = -sig*beta^(j)*s_d*d_d'*d_d        ; % Re-calculo la ganancia
%    x_d_next_j = x_d_act + beta^(j)*s_d*d_d ; % Siguiente posicion en la iteracion
%    f_d_next_j = (norm(A*x_d_next_j -b))^2  ; % Valor funcional del primer intento de armijo
%    gan_real = f_d_act - f_d_next_j         ;
%  end 
%  x_d_next  = x_d_next_j                    ; % Defino la posicion
%  grad_d_next = (2*(A*x_d_next - b)'*A)'    ; % Calculo el nuevo cetor de gradiente
%  % Almaceno valores  
%  x_d_hist = [x_d_hist, x_d_next]           ;
%  grad_d_hist = [grad_d_hist , grad_d_next] ;
%  % Calculo variacion
%  f_d_next = f_d_next_j                     ;
%  var_d = norm(f_d_next-f_d_act)/norm(f_0)  ;
%  var_d_hist = [var_d_hist, var_d]          ;
%  error_d_it = norm(x_min-x_d_next)/norm(x_min) ;
%  rel_error_d = [rel_error_d,error_d_it]        ;
%  iterations_d = [iterations_d,it_d]            ;
%  end
%
%figure
%hold on
%plot(iterations_c,rel_error_c,'*-')
%plot(iterations_d,rel_error_d,'*-')
%xl = xlabel('Numero de iteracion')
%yl = ylabel('Error relativo')
%lg = legend('Paso fijo','Paso exacto','Regla de armijo')
%grid minor
%set (xl, "fontsize", 16);
%set (yl, "fontsize", 16);
%set (lg, "fontsize", 16);
%
figure
plot(iterations_b,f_b_hist,'*-')
xl = xlabel('Numero de iteracion');
yl = ylabel('Valor funcional');
lg = legend('Paso decreciente');
grid minor
set (xl, "fontsize", 16);
set (lg, "fontsize", 16);
set (yl, "fontsize", 16);

x_fact = [0:-0.001:-R];
y_fact = -sqrt(R^2 .- x_fact.^2);

figure
tit = title('')
plot(x_b_hist(1,:), x_b_hist(2,:),'-*',-1/2,0,'r*',x_fact,y_fact,'k',x_fact,-y_fact,'k')
axis([-1/2 1/4 -1/4 1/2])
xl = xlabel('X');
yl = ylabel('Y');
lg = legend('Evolucion del algoritmo','Solucion en R^2','Limite del espacio factible');
grid minor
set (xl, "fontsize", 16);
set (yl, "fontsize", 16);
set (lg, "fontsize", 16);


figure
tit = title('')
plot(rel_error_b, 'b-*')
xl = xlabel('Iteracion');
yl = ylabel('Diferencia entre pasos consecutivos');
grid minor
set (xl, "fontsize", 16);
set (yl, "fontsize", 16);
