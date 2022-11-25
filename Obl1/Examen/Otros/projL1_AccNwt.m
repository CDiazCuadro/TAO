function[x, l, loops, dummy] = projL1AccNwt(b, tau, nMaxIter, flagR)
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

if nargin < 4
  flagR = 0;
  if nargin < 3
    nMaxIter = 20;
  end
end


% only for test_projL1.m purposes
dummy = 0;


if flagR == 0
  adaptR = nMaxIter + 1;
else
  adaptR = 2;
end

s0 = sign(b);
bnorm = s0'*b;
bAbs = s0.*b;


sN0 = length(b(:));


xnorm = bnorm;


if bnorm <= tau
  x = b;
  return;
end

% ==================================================
%  Init

      l0 = (bnorm - tau)/sN0;      
            
      s = s0.*(bAbs>l0);
      sN = s'*s;
      sb = s'*b;
      
      
      
      
% =================================

r = 1;

for k=1:nMaxIter

      
    if sN == sN0
      break;
    end    
    
    alpha = ( (sb - tau)/(sN) )/l0;

    if( k >= adaptR)
      r = 0.5/(alpha - 1.0);
      r = min( [r, 6] );
      r = max( [r, 1] );
    end
    
    l = l0*(1 - ((r-1)/r)*alpha) / ( (r+1)/r - alpha);
    
    sN0 = sN;
    
    s = s0.*(bAbs>l);
    sN = s'*s;
    
    if sN == 0
      l = ( (sb - tau)/(sN0) );
      s = s0.*(bAbs>l);
      sN = s'*s;      
    end
    
    
    sb = s'*b;
                             
    l0 = l;
      
end

l = ( (sb - tau)/(sN) );
    
x = s.*max(0, bAbs - l);  
loops = k;

    