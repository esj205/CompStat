gibbs_r <- function(N, thin){
  mat <- matrix(nrow=2, ncol=N)
  x <- y <- 0
  for (i in 1:N){
    for (j in 1:thin) {
      x <- rgamma(1, 3, y*y +4)
      y <- rnorm(1, 1 / (x+1), 1 / sqrt(2*(x+1)))
    }
    mat[,i] <- c(x,y)
  }
  mat
}

output = gibbs_r(10000, 5)
hist(output[1,], nclass=100)
hist(output[2,], nclass=100)


##Rcpp 이용하기

sourceCpp('gibbs_cpp.cpp')
output2 = gibbs_cpp(10000, 5)
hist(output2[1,], nclass=100)
hist(output2[2,], nclass=100)

sourceCpp('gibbs_cpp2.cpp')
output3 = gibbs_cpp2(10000, 5)
hist(output3$mat[1,], nclass=100)
hist(output3$mat[2,], nclass=100)

output4 = gibbs_cpp3(10000, 5)
hist(output4$x, nclass=100)
hist(output4$y, nclass=100)