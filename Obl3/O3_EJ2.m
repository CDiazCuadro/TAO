close all
clear all
cvx_setup
c1 = 1 ;
c2 = 4 ;
R  = 30;
d3 = 10;
limite = 50;
p1=zeros(1,limite);
p2=zeros(1,limite);
p3=zeros(1,limite);
g1=zeros(1,limite);
g2=zeros(1,limite);
t =zeros(1,limite);
lambda=zeros(1,limite);
costo=zeros(1,limite);
axis=[0:limite];

for d2=axis
%planteo del problema dentro del entorno cvx
cvx_begin 
cvx_precision high
variable x(6);  
dual variable l
minimize(x(4) + max(0,4*(x(5)-40)))
subject to
    x(1)+x(3)==x(4) 
    x(2)+x(3)+x(5)==d2
    l:x(1)-x(2)==d3
    -x(1)-x(2)+x(3)==0
    x(1)<=R
    -x(1)<=R
    max(0,4*(x(5)-40))<=x(6)
    x(4)>=0
    x(5)>=0
    x(6)>=0
cvx_end

p1(d2+1)=x(1);
p2(d2+1)=x(2);
p3(d2+1)=x(3);
g1(d2+1)=x(4);
g2(d2+1)=x(5);
t(d2+1) =x(6);
lambda(d2+1)=l;
costo(d2+1)=x(4) + max(0,4*(x(5)-40));
end

figure
plot(axis, g1,axis, g2,axis, p2)
legend('g1','g2','p1')
xlabel('d2')
ylabel('Potencia MW')

figure
plot(axis, lambda)
legend('Lambda')
xlabel('d2')
ylabel('Valor del multiplicador')

figure
plot(axis, t)
legend('Variable slack t')
xlabel('d2')
ylabel('Valor de slack')