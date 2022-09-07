close all
clear all


x = 0:0.001:0.42;
figure (1)
hold on
for c=0:0.1:0.6
  y = sqrt(10^(-c) + x.^2);
  plot(x,y, 'LineWidth',5*(c+0.4))
endfor
h = legend({'c=0.0','c=0.1','c=0.2','c=0.3','c=0.4','c=0.5','c=0.6'});      
set (h, "fontsize", 16);
y2=1-x.^2;
y3=2.*x;
y45=[0.5, 0.5, 1]; x45 = [0.25, 0,0] ;
ancho =3;
plot(x,y2,'k', 'LineWidth',ancho,x,y3,'k', 'LineWidth',ancho,x45,y45,'k', 'LineWidth',ancho)
axis([0 0.5 0.5 1.1])
xl=xlabel('x')
set (xl, "fontsize", 16);
yl=ylabel('y')
set (yl, "fontsize", 16);
tit = title('Dominio de f(x,y) - Curvas de nivel')
set (tit, "fontsize", 16);
grid minor


