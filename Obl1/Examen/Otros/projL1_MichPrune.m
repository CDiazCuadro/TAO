function[x, l, loops, dummy] = projL1_MichPrune(b, tau, nMaxIter)
%  
%  
% Author
% ======
% 
% Paul Rodriguez   prodrig@pucp.pe
%  
%  
% Related papers
% ==============
%  
%  [1] Paul Rodriguez, 'An accelerated Newton's method for projections onto the l1-Ball,' IEEE MLSP, 2017.
%  
%  [2] Paul Rodriguez, 'Accelerated Gradient Descent Method for Projections onto the l1-Ball', IEEE IVMSP 2018.
%  

if nargin < 3
  nMaxIter = 20;
end

l = 0;
loops = 0;

if tau==0
    x = b*0;
    return;
end

% only for test_projL1.m purposes
dummy = 0;


s0 = sign(b);
bnorm = s0'*b;
bAbs = s0.*b;

v = bAbs;

sN0 = length(b(:));


xnorm = bnorm;


if bnorm <= tau
  x = b;
  return;
end

% ==================================================
%  Init

      l0 = (bnorm - tau)/sN0;      
            
      v = v(v>l0);  

      sN = length(v);
      sb = sum(v);

      
      
      
      
% =================================

for k=1:nMaxIter

      
    if sN == sN0
      if(k==1) l = l0; end
      break;
    end    
    
    l = ( (sb - tau)/(sN) );
          
    sN0 = sN;
    
                             
    v = v(v>l);  
    sN = length(v);
    sb = sum(v);
      
end

    
x = s0.*max(0, bAbs - l);  
loops = k;
    