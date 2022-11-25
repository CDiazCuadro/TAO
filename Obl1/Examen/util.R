
proj_l1<-function(x){ 
  if(sum(abs(x)) < 1 ){
      return(x)
  }
  n <- length(x)
  a <- abs(x)
  s <- sign(x)
  idx <- order(a, decreasing = T)  
  y <- a[idx]
  w <- (cumsum(y) - 1 ) / seq(1,n)  
  j <- 0
  z <- 1
  while(z  > 1e-8){
    j <- j+1
    y = pmax(a - w[j],0)
    z = abs(sum(y) - 1)
    if( z < 1e-8){
      xp = s * y
      return(xp)
    }
  }
  print("no deberiamos llegar hasta aca")
}

x <- c(0.1,-0.7,0.8,-1)
y = proj_l1(x)
print(y)
