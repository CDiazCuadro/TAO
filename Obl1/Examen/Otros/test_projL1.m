function[lStats] = test_projL1(N, tau, loops, nMaxIter)
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

if nargin < 4
  nMaxIter = 20;
  if nargin < 3
    loops = [];
    if nargin < 2
      tau = 1;
      if nargin < 1
        N = 1e6;
      end
    end
  end
end


off = 1/N;

loops = setLoops(N, loops);

%  Some general definitions

pL1{1} = @(var1, var2) projL1_Condat(var1, var2);
pL1{2} = @(var1, var2) projL1_Duchi(var1, var2);
pL1{3} = @(var1, var2) projL1_Mich(var1, var2, nMaxIter);
pL1{4} = @(var1, var2) projL1_MichPrune(var1, var2, nMaxIter);
pL1{5} = @(var1, var2) projL1_AccNwt(var1, var2, nMaxIter, 0);
pL1{6} = @(var1, var2) projL1_AccNwt(var1, var2, nMaxIter, 1);
pL1{7} = @(var1, var2) projL1_Ntrv(var1, var2, nMaxIter);

L = length(pL1);


% ===================================
% ===================================

% experiment 1: mean = 1/N and SD = 1

a = 1;

for k = 1:loops

  b = a*randn(N,1) + off;

  for n=1:L
  
    t = tic;
    [x{n} l{n}(k) loops1{n}(k), loops2{n}(k)] = pL1{n}(b, tau);
    tCP{n}(k) = toc(t);
    l1Constrain{n}(k) = sum(abs(x{n})) - tau;
    
  end % _END_ FOR(n)
  

end % _END_ FOR(k)


lStats(1:L,:) = comp_dply_stats(L, 0, tCP, l, loops1, loops2, l1Constrain);


% ===================================
% ===================================

% experiment 2: mean = 1/N and SD = 1e-3

a = 1e-3;

for k = 1:loops

  b = a*randn(N,1) + off;

  for n=L+1:2*L
      
    t = tic;
    [x{n} l{n}(k) loops1{n}(k), loops2{n}(k)] = pL1{n-L}(b, tau);
    tCP{n}(k) = toc(t);
    l1Constrain{n}(k) = sum(abs(x{n})) - tau;
    
  end % _END_ FOR(n)
  

end % _END_ FOR(k)

lStats(L+1:2*L,:) = comp_dply_stats(L, 1, tCP, l, loops1, loops2, l1Constrain);



% ===================================
% ===================================

% experiment 3: mean = 0 and SD = 1e-3, except
%  one element at a random position, which is a 
%  random Gaussian number of mean 1 and SD 1e-3 .


a = 1e-3;

for k = 1:loops

  b = a*randn(N,1);
  
  bn = a*randn(1,1) + 1;
  n = randi([1,N],1);
  b(n) = bn;

  
  for n=2*L+1:3*L
  
    t = tic;
    [x{n} l{n}(k) loops1{n}(k), loops2{n}(k)] = pL1{n-2*L}(b, tau);
    tCP{n}(k) = toc(t);
    l1Constrain{n}(k) = sum(abs(x{n})) - tau;
    
  end % _END_ FOR(n)
  

end % _END_ FOR(k)

lStats(2*L+1:3*L,:) = comp_dply_stats(L, 2, tCP, l, loops1, loops2, l1Constrain);


return;


%  =======================================
%  =======================================

function[lStats] = comp_dply_stats(L, Loff, tCP, l, loops1, loops2, l1Constrain)

algos  = {'Condat    ', 'Duchi     ', 'Mich.     ', 'Mich+p    ', 'AccNwt    ', 'AccNwt+var', 'Nesterov  '};
Header = '              mean(T)    std(T)    lambda    std(lambda)   loops1   loops2   l1Cond (1e-9)';
sFMT = '%s  %9.2e %9.2e  %9.2e    %9.2e      %2.1f     %2.1f      %9.2e';
sFMT_bs = '%s  %9.2e %9.2e  %9.2e    %9.2e      \b%2.1f     %2.1f      %9.2e';

k = Loff*L;
% Compute stattictics
for n=1:L

  lStats(n,1) = mean(tCP{n+k});
  lStats(n,2) = std(tCP{n+k});
  lStats(n,3) = mean(l{n+k});
  lStats(n,4) = std(l{n+k});
  lStats(n,5) = mean(loops1{n+k});
  lStats(n,6) = mean(loops2{n+k});
  lStats(n,7) = std(l1Constrain{n+k});
  
end % _END_ FOR(k)


% Display stattictics
disp(' ');
disp(' ');
switch(Loff)
  case{0}
    disp('Case I');
  case{1}
    disp('Case II');
  case{2}
    disp('Case III');  
end
disp('-----------');
disp(' ');

disp(Header);
for n=1:L
 
  if lStats(n,5) < 10
    disp(sprintf(sFMT, ...
        algos{n}, lStats(n,1), lStats(n,2), lStats(n,3), lStats(n,4), ...
        lStats(n,5), lStats(n,6), (1e9)*lStats(n,7)));
  else
    disp(sprintf(sFMT_bs, ...
        algos{n}, lStats(n,1), lStats(n,2), lStats(n,3), lStats(n,4), ...
        lStats(n,5), lStats(n,6), (1e9)*lStats(n,7)));  
  end
  
  
end % _END_ FOR(n)

return;

%  ===========================================
%  ===========================================

function[loops] = setLoops(N, loops)

if isempty(loops)

  switch(N)
  
    case{ 1e5 }
      loops = 50;
    case{ 1e6 }
      loops = 25;
    case{ 1e7 }
      loops = 10;
    
    otherwise
      loops = 100*(N<1e5) + 20*(N>=1e5);
      
  end
end

return;


