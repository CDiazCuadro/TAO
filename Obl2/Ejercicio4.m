## Ejercicio 4

clear all, close all, clc
# Carga de datos
cd ./datos
load A.txt 
load B.txt
load C.txt 

load alpha.txt 
load beta.txt 
load gamma.txt

cd ..

## Defino matrices auxiliares D, delta y H
n = length(A);

# Matriz D
D = zeros (3*n,3*n);
D(1:n,1:n)              = A;
D(n+1:2*n,n+1:2*n)      = B;
D(2*n+1:3*n,2*n+1:3*n ) = C;

# Matriz H 
H = zeros(2*n,3*n);
H(1:n,1:n)            = eye(n);
H(1:n,n+1:2*n)        = -eye(n);
H(n+1:2*n, 1:n)       = eye(n);
H(n+1:2*n, 2*n+1:3*n) = -eye(n);

# Vector delta
delta = zeros(3*n,1);
delta = [alpha; beta; gamma];  

## Calculo solucion exacta x_op = inv(M)*m
# Defino M, m y mu
M  = A'*A + B'*B + C'*C;
m  = A'*alpha + B'*beta + C'*gamma;
mu = alpha'*alpha + beta'*beta + gamma'*gamma;
x_op = inv(M)*m

## Calculo loa multiplicadores lambda_op
lambda_op = inv(H*H')*H*D'*delta;

## Solucion Nummerica:
diff_max = 1e-8;
k_max = 100 ;
t0 = 1/norm(D);
w0 = rand(30,1);


# Parte e) lambda = lambda_op
#          tk = t02^k.

lambda = lambda_op ;
k      = 0          ;  
diff  = diff_max + 1;
w_hist  = [w0]      ; 
diff_hist = [diff]  ;

while diff > diff_max && k < k_max
      k = k + 1 ;
      t = t0*2^k;
      w_old = w_hist(:,end);
      w = inv( (1/2)*(D'*D + D*D') + t*H'*H)*(D'*delta - H'*lambda) ;
      diff = norm(w - w_old)/norm(w);
      diff_hist = [diff_hist,diff]  ;
      w_hist = [w_hist,w];
      
      
  
end
x_num = (w(1:10) + w(11:20) + w(21:30))/3