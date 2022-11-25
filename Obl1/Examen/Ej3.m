%%% EJERCICIO 3
%
close all, clear all
parte_a = 0;
parte_b = 1;
%
% ai) Continuidad de h(e)
% FUNCION h(e) con e mayores que 0
if parte_a
    e_p  = [0:0.1:50];
    h_p  = e_p.^2./(e_p+1);
    dh_p = (e_p.^2 + 2*e_p) ./(e_p+1).^2 ;
    %
    % FUNCION h(e) con e menores que 0
    e_n=[-50:0.01:0];
    h_n=e_n.^2./(-e_n+1);
    dh_n = (3*e_n.^2+2*e_n) ./ (-e_n+1).^2 ;
    %
    % FUNCION error cuadratico error2(e) = e^2
    e_t = [e_n,e_p];
    error2 = e_t.^2;

    figure
    plot(e_n, h_n,e_p, h_p, e_t, error2)
    legend('h(e) con e<0','h(e) con e>0','e^2')
    grid minor
    xlabel('Valores de e entre [-1 1]')
    ylabel('h(e) - e^2')
    axis([-1 1 0 1])

    figure
    plot(e_n, h_n,e_p, h_p, e_t, error2)
    legend('h(e) con e<0','h(e) con e>0','e^2')
    grid minor
    xlabel('Valores de e entre [-50 50]')
    ylabel('h(e) - e^2')

    figure (2)
    plot(e_n, dh_n, e_p, dh_p)
end

%%%%%%%%%%%%%%%%%%%%%%%
if parte_b
    % Loading data
    load A_May.txt 
    A = A_May;
    load b.txt
    load x_robusto.txt
    x_op = x_robusto;
    %
    % -> Descenso por gradiente
    %     > x_next = x_act - beta^j.gradF(x_act)
    %     > Siendo gradF = A * grad_h(e)
    %
    % --> d) regla de Armijo    
    %      - Defino: sig = 0.1
    %                beta = 1/5       
    %      --> iterar j hasta que: f_old - f_next > - sig*b^j*s*gradF*dk
    %
    %
    % Parámetros de control del programa:
    tol_f = 1e-5                ; % Tolerancia para f(x)
    var = tol_f + 1             ; % Aseguro que var > tol la primera vez
    it_max = 5000               ; % Determino el maximo de iteraciones
    sig = 0.1                   ; % Parametro 1 de armijo
    beta = 1/5                  ; % Parametro 2 de armijo
    s = 1                       ; % Paso original
    % Variables iniciales:
    n = size(A,2)               ;
    rand ("seed", pi())         ; % Para que x0 sea siempre la misma
    x0 = 1*rand(n,1)          ; % Elijo un primer punto de Rn
    vec0 = A*x0-b               ; % Defino el vector a evaluar en h inicial
    f0  = sum(vec_h(vec0))      ; % Función de costo evaluado en x0
    var_0 = f0                  ; % Defino la cariacion de la f como la f inicial
    grad0 = A*grad_h(vec0)     ; % Gradiente inicial de la funcion de costo
    it = 0                      ; % Inicializo las iteraciones
    % Inicializo vectores de guardado de informacion
    X = [x0]                    ; % inicializo la matriz de historial de posiciones
    F = [f0]                    ; % Inicializo el vector historial de la funcion de costo 
    GRAD = [grad0]              ; % inicializo la matriz de almacenamiento de gradientes
    VAR = [var_0]               ; % inicializo el vector de historial de variacion
    real_error_0 = norm(x_op-x0); % Calculo el error inicial
    ERROR = [real_error_0]      ; % inicializo el vector de historial de errores relativos
    IT = [0]                    ;


    while var > tol_f && it < it_max
      it = it + 1                           ; % Actualizo el numero de iteracion
      d       = -GRAD(:,end)                ; % Defino direccion opuesto al gradiente
      x_act   = X(:,end)                    ; % Obtengo la posicion actual
      vec_act = A*x_act-b                   ; % Calculo el vector a evaluar h
      f_act = sum(vec_h(vec_act))           ; % Valor funcional de la posicion actual
      j = 1                                 ; % iterador de armijo
      gan = -sig*beta^(j)*s*d'*d            ; % Ganancia si fuera lineal
      x_next_j = x_act + beta^(j)*s*d       ; % Siguiente posicion en la iteracion j=1
      vec_next_j = A*x_next_j-b             ; % Calculo el vector a evaluar h
      f_next_j = sum(vec_h(vec_next_j))     ; % Valor funcional del primer intento de armijo
      gan_real = f_act - f_next_j           ;
      while  gan_real < gan
        j = j + 1                           ; % Aumento el iterador de armijo
        gan = -sig*beta^(j)*s*d'*d          ; % Re-calculo la ganancia
        x_next_j = x_act + beta^(j)*s*d     ; % Siguiente posicion en la iteracion
        vec_next_j = A*x_next_j-b           ; % Calculo el vector a evaluar h
        f_next_j = sum(vec_h(vec_next_j))   ; % Valor funcional del primer intento de armijo
        gan_real = f_act - f_next_j         ;
      end  
      x_next  = x_next_j                    ; % Defino la ultima it como la posicion siguiente definitiva
      error_real = norm(x_op-x_next)        ; % Calculo el error real con respecto al optimo
      vec_next = vec_next_j                 ; % Declaro el vector a evaluar h siguiente
      f_next = f_next_j                     ; % Defino con la ultima it el valor funcional de x_next
      grad_next = A*grad_h(vec_next)       ; % Calculo el nuevo vector de gradiente
      var = norm(f_next-f_act)              ;
      % Almaceno valores  
      X = [X, x_next]                       ;
      GRAD = [GRAD , grad_next]             ;
      F = [F, f_next]                       ;
      VAR = [VAR, var]                      ;
      ERROR = [ERROR,error_real]            ;
      IT = [IT,it]                          ;
    end
end