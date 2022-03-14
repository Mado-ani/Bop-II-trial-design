

 #### I made some adjustment on my code to estimate the sample size
 #### under the null and alternative hypothesis in same way
 #### and present the mean and standard deviation for both 
 #### values of theta as a matrix.

set.seed(201408673)
library(OpenMx)

lambda <- 0.5

gamma <- 0.8 

n1 <- 35  

n2 <- 55

# set the ( seed) number to get the same results when running the code agin.

# Set the function for our evaluation design.

ptm <- proc.time()

Bop2_estimation = function(lambda, gamma, n1, n2){
  
  Total = 10^4 
  
  theta <- c(0.5, 0.7)
  
  sample_simulate <- matrix(NA, Total,2)
  
  G <- matrix(NA, nrow = length(theta), ncol = 2)
  
  Response = rep(NA, length(theta))
  
  P_futility = rep(NA, length(theta))
  
  for (i in 1:Total) {
    for (j in 1:length(theta)) {
      
      
      
      
      Response[j] <- rbinom(1, n1, theta[j])
      
      
      a1 <- 0.5 + Response[j]
      
      b1 <- 0.5 + n1 - Response[j]
      
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

 #### in this code we  will obtain the estimated sample size
 #### under both the null and alternative hypotheses.
 ### the results will appear as a matrix 

    Bop2_estimation(lambda, gamma, n1, n2)  

   #                       [mean]      [SD]
   # [H_0 : theta = 0.5]   47.616   0.09652246
   # [H_1 : theta = 0.7]   54.946   0.01037878
    
    
 ##################
 ## counting the required time to perform the code of BOP2 design.

    ptm <- proc.time() ### starting the counter 

    y <- replicate(10, Bop2_estimation(lambda, gamma, n1, n2))  ### repeating the code for 10 times.

    proc.time() - ptm   ### end the counter.

 ### the results of our counter.
 ###  user  system elapsed 
 #    0.95    0.00    0.95 

############################################
    
    ## this is the code of estimating the sample size by using Monte Carlo Method
    ## without any change of the original code from the lecture notes..
    
    BOP2_estimation = function(lambda, gamma, n1, n2, theta)
     {
      Total = 10^4 
      
      sample_simulate <- rep(NA, Total)
      
      for (i in 1:Total) 
        {
         Responce_1 <- rbinom(1, n1, theta)
        
         a1 <- 0.5 + Responce_1
        
         b1 <- 0.5 + n1 - Responce_1
        
         P_futility <- pbeta(0.5, a1, b1)
        
         Cutoff <- 1 - lambda * (n1 / n2) ^ gamma
         
         if (P_futility > Cutoff) 
           {
            sample_simulate[i] <- n1
           } 
         else
           {
            sample_simulate[i] <- n2
           }
        }
      
      return(mean(sample_simulate))
     }    
     
    
    BOP2_estimation(lambda, gamma, n1, n2,0.5) 
    # [H_0 : theta = 0.5] [47.616]
    
    BOP2_estimation(lambda, gamma, n1, n2,0.7) 
    # [H_0 : theta = 0.7] [54.946]
    
    
    ##############
    ## calculating the time for running this code
    ptm <- proc.time() ### starting the counter
    
    
    Y= replicate (1,BOP2_estimation(lambda, gamma, n1, n2,0.5))
    
    BOP2_estimation(lambda, gamma, n1, n2,0.7)
    
    proc.time() - ptm   ### end the counter.
    
    ### this is the whole time of executing the BOP2_estimation code 
    ### for one time
    # user   system  elapsed 
    # 0.14    0.00     0.14
                     
    
    ptm <- proc.time() ### starting the counter
    
    
    Y= replicate (20,BOP2_estimation(lambda, gamma, n1, n2,0.5))
    
    BOP2_estimation(lambda, gamma, n1, n2,0.7)
    
    proc.time() - ptm   ### end the counter.
    
    ## while when we executed it for 20 times it took
    # user   system  elapsed 
    # 1.34    0.00     1.34
    
    ####################################################
    
    ## trying to improving the efficacy of our code 
    ## this code of finding the effective sample size.
    
    # finding the probability of responses in stage 1 conditional to theta
    # with first group of patients n1.
    
    prob_Response <- function(Response, n1, theta)
     {
       dbinom(Response,n1,theta)
     }
    
    BOP2_improve_eff <- function(lambda, gamma, n1, n2,theta) 
    {
      
      Response1_s<- 0: n1
      
      Cutoff1<- 1 - lambda * (n1 / n2)^gamma
      
      quits<- pbeta(0.5, Response1_s + 0.5,  n1 - Response1_s + 0.5)  > Cutoff1
      
      Probability_Response<- prob_Response(Response1_s, n1, theta)
      
      sum(n1 *quits *  Probability_Response  +  n2 * (!quits) * Probability_Response)
    }
      
    
    BOP2_improve_eff(lambda, gamma, n1, n2,0.5)
    
    # [H_0 : theta = 0.5]  [47.64121]
    
    BOP2_improve_eff(lambda, gamma, n1, n2,0.7)
    
    # [H_1 : theta = 0.7]  [54.95342]
    
    
    ## the results that appeared are approximately equal to the results
    ## that obtained from Monte Carlo estimation.
    
  ############################################################  
    
    ptm <- proc.time() ### starting the counter
    
    
    Y= replicate (20,BOP2_improve_eff(lambda, gamma, n1, n2,0.5))
    
    Y= replicate (20,BOP2_improve_eff(lambda, gamma, n1, n2,0.7))
    
    proc.time() - ptm   ### end the counter.
    
    ## this result presents that the improved code {BOP2_improve_eff}
    ## spend time less than the estimation code {BOP2_estimation}
    ## to estimate the expected sample size.
    
    # user  system elapsed 
    # 0.00    0.02    0.02
    
 ################################################################
    
   ## Finding the optimal value for lambda and gamma
   ## we set the counter of the values of lambda and gamma will change at each one decimal.
   ## because when we implemented it for more than that it took more time.

    ptm <- proc.time() ### starting the counter
    
    eval <- expand.grid(lambda = seq(0, 1, 0.1),gamma = seq(0, 1, 0.1))
      
    expected_sample_size <- rep(NA, nrow(eval))
      
        for(i in 1: nrow(eval))
          {
           expected_sample_size[i] <- Bop2_estimation(eval[i, 1], eval[i, 2], n1, n2)[1]
          }
    proc.time() - ptm   ### end the counter.

   ### the results of our counter.   
   ##    user  system elapsed 
   ##    11.55    0.02   11.55 
  
    
    expected_sample_size 
    
   ##   the minimum expected value is  
    min(expected_sample_size)

   # [35]
    
  #############################################################
   ### Testing the code that the sum of probabilities should equal to 1.
   ### by calculating the probability of observing the responses 

    probability_responses= function(Responses, n1)
     {
      choose(n1 , Responses) * beta(Responses+0.5, n1 - Responses + 0.5) / beta(0.5, 0.5)
     }
    probability_responses(10, 35)
    
    
    Test_probability_responses = function()
     {
      n1 = 35
      tot = 0
      for (z in 0:n1) 
        { 
          tot = tot + probability_responses(z, n1)
        
        }
      return( tot )
     }
 
    Test_probability_responses()
  
    ### the sum of probabilities is equal to 1.
    ### [tot] 1

   #######################################
    
    #### estimate the sample sizes under different values of theta 
    ### to investigate it is sensitivity 
    
    sensivity = function(lambda, gamma)
      {
       values_theta<- expand.grid(theta = seq(0.1, 0.9, 0.1))	
      
       Snew_expected <- rep(NULL, nrow(values_theta))	  
      
       for(e in 1: nrow(values_theta)) 
         {
        
         Snew_expected[e] <- BOP2_improve_eff(lambda, gamma, n1 = 30, n2 = 70, values_theta[e,1])
         }
      
       return(cbind(values_theta ,Snew_expected))
      
     }
    sensivity(lambda, gamma)
    
    ###  theta    Snew_expected
    # 1   0.1       30.00001
    # 2   0.2       30.03607
    # 3   0.3       31.60210
    # 4   0.4       41.41982
    # 5   0.5       58.30671
    # 6   0.6       68.07553
    # 7   0.7       69.91501
    # 8   0.8       69.99958
    # 9   0.9       70.00000
    
    
    
    #### Plotting the expected sample sizes ageist 
    #### the different values of null hypothesis and with different values 
    #### of the parameters gamma and lambda.  
    #### 

    par(mar=c(2, 1, 1, 1))
    ## Set a sequence for the values of theta.
    seq_theta <-seq(0.1, 0.9, 0.1) 
    
    ## To plot the expected samples sizes Vs the null hypothesis
    plot(seq_theta, sensivity(0.9, 0.2), type ="b", ylab= "The expected Sample Size",
         xlab=" The Null Hypothesis", lwd=3,lty = 2, col="blue") 
    
    ## Adding some lines for different values of gamma and lambda.  
    lines(seq_theta, sensivity(0.4, 0.9), type ="b", lwd=3,lty = 2, col="orange")  
    lines(seq_theta, sensivity(0.5, 0.2), type ="b", lwd=3,lty = 2, col="skyblue") 
    
    legend(.63, 38, c('(lambda = 0.9 , gamma = 0.2)' ,
                      '(lambda = 0.4 , gamma = 0.9)' ,
                      '(lambda = 0.5 , gamma = 0.2)'),
           lty=c(2,2,2), col=c('blue'  ,
                               'orange',
                               'skyblue'))
    





