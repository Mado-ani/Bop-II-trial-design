




G <- matrix(NA, nrow = n, ncol = n)

# Fill in the elements of G:

for(i in 1:n){
  for(j in 1:n){
    
    G[i, j] <-  as.numeric(exam$Group[i] == exam$Group[j]) 
  }  
}

set.seed(201408673)
library(OpenMx)

lambda <- 0.5

gamma <- 0.8 

n1 <- 30  

n2 <- 60

# set the ( seed) number to get the same results when running the code agin.





# Set the function for our evaluation design.

Bop2_evaluation = function(lampda, gamma, n1, n2){
  
  
  Total = 10^4 
  
  
  theta <- c(0.5, 0.7)
  
  
  sample_simulate <- matrix(NA, Total,2)
  
  G <- matrix(NA, nrow = length(theta), ncol = 2)
  Responce = rep(NA, length(theta))
  P_futility = rep(NA, length(theta)) 
  for (i in 1:Total) {
    for (j in 1:length(theta)) {
      
      
      
      
      Responce[j] <- rbinom(1, n1, theta[j])
      
      
      a1 <- 0.5 + Responce[j]
      
      b1 <- 0.5 + n1 - Responce[j]
      
      P_futility[j] <- pbeta(0.5, a1, b1)
      
      Cutoff <- 1 - lambda * (n1 / n2) ^ gamma
      
      if (P_futility[j] > Cutoff) {
        
        sample_simulate[i,j] <- n1
        
      } else {
        
        sample_simulate[i,j] <- n2
        
      }
    }
  }
  return(matrix(c(colMeans(sample_simulate), t(diag2vec(sqrt(var(sample_simulate )/ Total* diag(2) )))),2,2))
}    

Bop2_evaluation(lampda, gamma, n1, n2)  
