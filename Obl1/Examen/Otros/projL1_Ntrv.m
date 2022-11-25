function[x, l, loops, loops2] = projL1_Nvr(b, tau, nMaxIter)
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
%  [1] Paul Rodriguez, 'Accelerated Gradient Descent Method for Projections onto the l1-Ball', IEEE IVMSP 2018.
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


s0 = sign(b);
bAbs = abs(b);

bnorm = sum(bAbs);

v = bAbs;

N0 = length(b(:));


if bnorm <= tau
  x = b;
  return;
end

% ==================================================
%  Init

t0 = 1;
t1 = 0.5*(1 + sqrt(1 + 4*t0*t0));

l0 = ( bnorm-tau )/N0;

v = v(v>l0);
N = length(v);

sb = sum(v);

% ------------------

for k = 1:nMaxIter,

  lMich = (sb-tau)/N;
  lNwt  = l0/(2.0 - lMich/l0);

  
  t0 = t1;
  t1 = 0.5*(1 + sqrt(1 + 4*t0*t0));
  lNrv = lNwt + ((t0-1)/t1)*(lNwt - lMich);

  
  
  N0 = N;
  v = v(v>lMich);
  zz = v(v>lNrv);
  sb = sum( zz );
  N = length(zz);
  
  if( -sb + lNrv*N + tau > 0 )
    l = lMich;
    break;
  else
    l = lNrv;
    v = zz;
  end
  
  if k==1
    l0 = lNwt;
  else
    l0 = lMich;
  end
  
  
end

loops = k;

%  =======================


N = length(v);
sb = sum(v);

N0 = length(b);

for k = 1:nMaxIter,

    if N == N0
      break;
    end    
    
    l = ( (sb - tau)/(N) );
          
    N0 = N;
    
                             
    v = v(v>l);  
    N = length(v);
    sb = sum(v);

end

loops2 = k;

%  =======================

x = s0.*max(0, bAbs - l);  

