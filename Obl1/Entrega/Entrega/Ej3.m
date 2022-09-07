#### Ejercicio 3
close all
clear all
clc

%% a)
x = [-1:0.01:1];
y = 4*x.^4 - x.^3 -4*x.^2 +1;

func = [x,y]';
[y_min,iy_min] = min(y);
x_min = x(iy_min); 

figure (1)
plot(x,y, x_min,y_min,'o', 'LineWidth',3) 
grid minor
tit = title('EJ 3a) - funcion f(x,y) y minimo global');
set (tit, "fontsize", 16);
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
h = legend({'f(x)','Mínimo'});      
set (h, "fontsize", 16);

%% b)
x = [-1:0.01:1];
y = x.^3;

func = [x,y]';
[y_min,iy_min] = min(y);
x_min = x(iy_min); 

figure (2)
plot(x,y, x_min,y_min,'o', 'LineWidth',3) 
grid minor
tit = title('EJ 3b) - Funcion f(x,y) y mínimo');
set (tit, "fontsize", 16);
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
h = legend({'f(x)','Mínimo'});      
set (h, "fontsize", 16);

%% c)

% Si a pertence a el dominio --> el minimo es en a
% Si a <=-1   --> el minimo es x=-1
% Si a >= 1  --> el minimo es x=1

x = [-1:0.01:1];
a = -0.5;
y = (x-a).^2 + 1;

func = [x,y]';
[y_min,iy_min] = min(y);
x_min = x(iy_min); 

figure (3)
plot(x,y, x_min,y_min,'o', 'LineWidth',3) 
grid minor
tit = title('EJ 3c) - funcion f(x,y) y mínimo');
set (tit, "fontsize", 16);
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
h = legend({'f(x)','Mínimo'});      
set (h, "fontsize", 16);


%% d)

% f(x1,x2): R2 --> R

xm = [2,1/2];
scale=1e1;
##x2a = [0:1/scale:1]';
##n=length(x2a);
##figure (4)
##hold on
##for x1a = 0:1/scale:1
##  x1_plot=x1a*ones(n,1);
##  for k = 1:n
##    xda=[x1a,x2a(k)];   
##    yda(k,1) = norm(xda-xm) + 1;
##  end
##  plot3(x1_plot,x2a,yda(:,1),'LineWidth',2*(2-x1a))
##end  
##
##x1b = x2a;
##for x2b = 0:1/scale:1
##  x2_plot=x2b*ones(n,1);
##  for k = 1:n
##    xdb=[x1b(k),x2b];   
##    ydb(k,1) = norm(xdb-xm) + 1;
##  end
##  plot3(x1b,x2_plot,ydb(:,1),'LineWidth',2*(2-x2b))
##end  

tx = ty= [0:1/scale:1];
tm= xm;
[xx, yy] = meshgrid (tx, ty);
for i = 1:length(tx)
  for j = 1:length(ty)
    txy = [tx(i), ty(j)];
    tz(i,j) = (norm(txy - tm))^2;
  end
end
figure (5)
mesh (tx, ty, tz)
hold
[min_y, i_min_y] = min(tz);
[min_t, i_min_x] = min(min_y);
plot3(tx(i_min_x),ty(i_min_y(1)),min_t,'o', 'LineWidth',3)
grid minor
tit = title('EJ 3d) Funcion f(x,y)');
set (tit, "fontsize", 16);
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
zl=zlabel('z');
set (zl, "fontsize", 16);
h = legend({'f(x)','Mínimo'});      
set (h, "fontsize", 16);

figure (6)
hold on
x2d = ty;
for c = -0.4:0.05:0
  x1 = 2 - sqrt( (c-1)^4 .- (x2d .- 1/2).^2 );
  plot(x2d,x1,'LineWidth',5*(0.5-c))
end
axis([0 1 0 1])
grid minor
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
tit = title('EJ 3d) Curvas de nivel en el dominio de f(x,y)');
set (tit, "fontsize", 16);


figure (7)
hold on
s_x = [-1,2];
f_s = [0:0.2:1];
for it = f_s
  plot(s_x,[it,it],'b','LineWidth',3)
end

for it = f_s
 plot([it,it],s_x,'r','LineWidth',3)
end
lg = legend('S_2','','','','','','S_1','','','','','')
grid minor
xl=xlabel('x');
set (xl, "fontsize", 16);
yl=ylabel('y');
set (yl, "fontsize", 16);
tit = title('EJ 3d) - intersección de S_1 y S_2');
set (tit, "fontsize", 16);
