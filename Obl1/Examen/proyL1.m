function xp = proyL1(x)

  if sum(abs(x)) < 1
      xp = x ;
  
  else
    n = length(x)               ;
    a = abs(x)                  ;
    s = sign(x)                 ;
    [y,idx] = sort(a,'descend') ; 
    seq = [1:n]                 ; 
    seq = reshape(seq,size(y))  ;
    w = (cumsum(y) - 1 ) ./ seq ;
    j = 0                       ;
    z = 1                       ;
    
    while z  > 1e-8
        j = j+1                     ;
        y = max([a - w(j);zeros(size(a))]) ;
        z = abs(sum(y) - 1)         ;
        if z < 1e-8
            xp = s .* y;              
        end %if
    end %while
  end
end

