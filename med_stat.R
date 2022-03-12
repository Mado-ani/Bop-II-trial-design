

   # In this code we try to simulate the Bayesian optimal phase II trial 
    # by implementing this method on the RESERVE clinical trial.....
                
                     ######################
                
    # Defining the enetial value for the parameters.
                     
    
     
    # Total           : is the total number that will used in the desigen.
    # [n1, n2]        : the number of participants that recruited in each stage( stage 1 , stage 2)
    # theta           : is a vector for the values of theta under the Null and Alternative hypothesis.
    # sample_simulate : is an empty vector to store simulated sample size .
    # [a1 , b1]       : the parameters for the posterior distribution Beta ~ (a1, b1)
    # Responce_1      : is the number of responce at stage one from the RESERVE trial. 
    # P_futility      : is the probaability of futility.
    # Cutoff          : is the threshold to determine progression, based on the decision rule.

# set the ( seed) number to get the same results when running the code agin.
set.seed(201408673)
library(OpenMx)

lambda <- 0.5

gamma <- 0.8 

n1 <- 30  

n2 <- 60


Bop2_evaluation = function(lampda, gamma, n1, n2){
  
  Total = 10^4 
  
  theta <- c(0.5, 0.7)
  
  sample_simulate <- matrix(NA, Total,2)      # define an empty matrix to save the results in.
  
  Responce = rep(NA, length(theta))           # define an empty vector to saved the responces values in.
  
  P_futility = rep(NA, length(theta))         # define an empty vector to saved the probability of futility in.
  
  for (i in 1:Total)
  {
    for (j in 1:length(theta))
    {
      Responce[j] <- rbinom(1, n1, theta[j])
      
      a1 <- 0.5 + Responce[j]
      
      b1 <- 0.5 + n1 - Responce[j]
      
      P_futility[j] <- pbeta(0.5, a1, b1)
      
      Cutoff <- 1 - lambda * (n1 / n2) ^ gamma
      
      if (P_futility[j] > Cutoff)
          {
           sample_simulate[i,j] <- n1
          } 
      else 
         {
          sample_simulate[i,j] <- n2
         }
    }
  }
  
  # this code will return the matrix that contain the mean and standard deviation for the expected values under the null (0.5) and the alternative hypothesis.
  return(matrix(c(colMeans(sample_simulate), t(diag2vec(sqrt(var(sample_simulate )/ Total* diag(2) )))),2,2))
}    

Bop2_evaluation(lampda, gamma, n1, n2)  
