function[x, l, dummy1 dummy2] = projL1Condat(b,c)
% Projection onto the simplex or the l1 ball.
%  
% Author: Laurent Condat, PhD, CNRS research fellow in Grenoble, France.
% Version: 1.0, Sept. 18, 2015.
% Copyright: Laurent Condat.
%  
%  
% Notes
% ======
%  
% Downloaded from https://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/download/proj_simplex_l1ball.m
%  
% Minimum modifications for its usage in
%  
%  [1] Paul Rodriguez, 'An accelerated Newton's method for projections onto the l1-Ball,' IEEE MLSP, 2017.
%  
%  [2] Paul Rodriguez, 'Accelerated Gradient Descent Method for Projections onto the l1-Ball', IEEE IVMSP 2018.
%  


if norm(b,1) <= c
  x = b;
  l = 0;
  return;
end

% only for test_projL1.m purposes
dummy1 = 0;
dummy2 = 0;



l = max(max((cumsum(sort(abs(b),1,'descend'),1)-c)./(1:size(b,1))'),0);

x = max(abs(b)-l,0).*sign(b);