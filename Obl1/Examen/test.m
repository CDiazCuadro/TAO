clear all, close all
xx= [-1:0.01:1];
yy = [0:0.01:2];
A=[2 1; -1 0];
b=[2;3];
i=0;
j=i;
XX=xx'*ones(size(xx));
YY=yy'*ones(size(yy));
for ix = xx
    i = i +1;
    j = 0;
    for iy = yy
        j = j +1;
        F(i,j) = (1/2)*norm(A*[ix;iy] - b)^2 + norm([ix;iy],Inf);
    end
end
plot3(XX,YY',F)
xlabel('X')
