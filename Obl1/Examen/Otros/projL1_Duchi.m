function[x, l, dummy1, dummy2] = projL1_Duchi(b, c)
% PROJECTONTOL1BALL Projects point onto L1 ball of specified radius.
%
% w = ProjectOntoL1Ball(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w||_1 <= b.
%
%   That is, performs Euclidean projection of v to the 1-norm ball of radius
%   b.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)
%  
%  
% Notes
% ======
%  
% Downloaded from  https://stanford.edu/~jduchi/projects/DuchiShSiCh08/ProjectOntoSimplex.m
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


u = sort(abs(b), 'descend');
sv = cumsum(u);

rho = find( u > (sv - c)./(1:length(u))', 1, 'last');
l = max(0, (sv(rho) - c) / rho);


x = shrink(b, l);


