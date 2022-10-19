close all
clear all
cvx_setup
c1 = 1 ;
c2 = 4 ;
R  = 30;
d3 = 10;

c=[0 0 0 c1 c2 0]';
A = [1 0 0 0 0   0;...
    -1 0 0 0 0   0;...
     0 0 0 0 c2 -1;...
     0 0 0 -1 0  0;...
     0 0 0  0 -1 0];
b=[R R 160 0 0]';
C = [1  0 1 -1 0 0;...
     0  1 1  0 1 0;...
     1 -1 0  0 0 0;...
    -1 -1 1 0 0 0];




limite = 50;
p1=zeros(1,limite);
p2=zeros(1,limite);
p3=zeros(1,limite);
g1=zeros(1,limite);
g2=zeros(1,limite);
t =zeros(1,limite);
lambda=zeros(4,limite);
mu=zeros(5,limite);
costo=zeros(1,limite);
axis=[0:limite];

for d2=axis
    d = [0 d2 d3 0]';
%planteo del problema dentro del entorno cvx
cvx_begin 
cvx_precision high
variable x(6);  
dual variable l
dual variable m
minimize(x(4) + max(0,4*(x(5)-40)))
subject to
    m:A*x <= b
    l:C*x == d
cvx_end

p1(d2+1)=x(1);
p2(d2+1)=x(2);
p3(d2+1)=x(3);
g1(d2+1)=x(4);
g2(d2+1)=x(5);
t(d2+1) =x(6);
lambda(:,d2+1)=l;
mu(:,d2+1) = m;
costo(d2+1)=x(4) + max(0,4*(x(5)-40));
end

figure
plot(axis, g1,axis, g2,axis, p2)
legend('g1','g2','p1')
xlabel('d2')
ylabel('Potencia MW')

figure
plot(axis, lambda(3,:))
legend('Lambda')
xlabel('d2')
ylabel('Valor del multiplicador')

figure
plot(axis, t)
legend('Variable slack t')
xlabel('d2')
ylabel('Valor de slack')