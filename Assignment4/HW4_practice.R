f = function(x, mean, sd){
  out = dnorm(x, mean, sd)*dcauchy(x, 5,2)
  return(out)
}

#riemann
riemann = function(interval, n, data, sd){
  h   = (interval[2]-interval[1])/n
  x   = interval[1] + (0:(n-1))*h
  out = h*sum(f(x, mean(data), sd))
  return(out)
}


1/riemann(interval, n=100, data=x, sd=sd)


#trapezoidal
trapezoidal = function(interval, n, data, sd){
  h   = (interval[2]-interval[1])/n
  x   = interval[1] + (1:(n-1))*h
  low = h/2*f(interval[1], mean(data), sd)
  middle = h*sum(f(x,mean(data),sd))
  high = h/2*f(interval[2], mean(data), sd)
  out = sum(low, middle, high)
  return(out)
}

1/trapezoidal(interval, n=100, data=x, sd=sd)

#simpson
simpsons = function(interval, n, data, sd){
  h   = (interval[2]-interval[1])/n
  x   = interval[1] + (0:n)*h
  out = NULL
  for(i in 1:(n/2)){
    out[i] = h/3*f(x[2*i - 1], mean(data), sd)
    out[i] = out[i] + 4*h/3*f(x[2*i],mean(data), sd)
    out[i] = out[i] + h/3*f(x[2*i + 1],mean(data), sd)
  }
  res=sum(out)
  return(res)
}

1/simpsons(interval, n=100, data=x, sd=sd)

f = function(x, mean, sd){
  out = dnorm(x, mean, sd)*dcauchy(x, 5,2)
  res = 7.84654*out
  return(res)
}


interval <- c(2,8)


riemann_table <- function(interval, x, sd, tol){
  estimation = c()
  relative_error = c()
  sub_interval = c()
  iter = 1
  n = 2
  diff = 1
  while(diff > tol){
    if(iter==1){
      first_est = riemann(interval, n, x, sd)
      estimation[iter] = first_est
      sub_interval[iter] = n
      relative_error[iter] = 0
      n=n+2
      iter = iter+1
      old_estimation = first_est
    }
    else{
      new_estimation=riemann(interval, n, x, sd)
      diff = abs(new_estimation - old_estimation)
      estimation[iter] = new_estimation
      sub_interval[iter] = n
      relative_error[iter] = new_estimation - old_estimation
      n=n+2
      iter = iter+1
      old_estimation = new_estimation
    }
  }
  res = data.frame(estimation=estimation, sub_interval=sub_interval, relative_error=relative_error)
  return(res)
}

riemann_table(interval, x, sd, tol=0.0001)


trapezoidal_table <- function(interval, x, sd, tol){
  estimation = c()
  relative_error = c()
  sub_interval = c()
  iter = 1
  n = 2
  diff = 1
  while(diff > tol){
    if(iter==1){
      first_est = trapezoidal(interval, n, x, sd)
      estimation[iter] = first_est
      sub_interval[iter] = n
      relative_error[iter] = 0
      n=n+2
      iter = iter+1
      old_estimation = first_est
    }
    else{
      new_estimation=trapezoidal(interval, n, x, sd)
      diff = abs(new_estimation - old_estimation)
      estimation[iter] = new_estimation
      sub_interval[iter] = n
      relative_error[iter] = new_estimation - old_estimation
      n=n+2
      iter = iter+1
      old_estimation = new_estimation
    }
  }
  res = data.frame(estimation=estimation, sub_interval=sub_interval, relative_error=relative_error)
  return(res)
}

trapezoidal_table(interval, x, sd, tol=0.0001)


simpsons_table <- function(interval, x, sd, tol){
  estimation = c()
  relative_error = c()
  sub_interval = c()
  iter = 1
  n = 2
  diff = 1
  while(diff > tol){
    if(iter==1){
      first_est = simpsons(interval, n, x, sd)
      estimation[iter] = first_est
      sub_interval[iter] = n
      relative_error[iter] = 0
      n=n+2
      iter = iter+1
      old_estimation = first_est
    }
    else{
      new_estimation=simpsons(interval, n, x, sd)
      diff = abs(new_estimation - old_estimation)
      estimation[iter] = new_estimation
      sub_interval[iter] = n
      relative_error[iter] = new_estimation - old_estimation
      n=n+2
      iter = iter+1
      old_estimation = new_estimation
    }
  }
  res = data.frame(estimation=estimation, sub_interval=sub_interval, relative_error=relative_error)
  return(res)
}


simpsons_table(interval, x, sd, tol=0.0001)


riem <- riemann_table(interval, x, sd, tol=0.0001)
trp <- trapezoidal_table(interval, x, sd, tol=0.0001)
simps <- simpsons_table(interval, x, sd, tol=0.0001)

riem
trp
simps

abs(riem[13,1]-0.99605)
