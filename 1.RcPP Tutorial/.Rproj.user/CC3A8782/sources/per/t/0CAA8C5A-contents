sourceCpp("example5.cpp")

output = gibbs_cpp(10000, 5)
hist(output$mat[,1],nclass=100)
hist(output$mat[,2],nclass=100)

output2 = gibbs_cpp2(10000, 5)
hist(output2$x,nclass=100)
hist(output2$y,nclass=100)
