library(xptr)
library(pracma)
source("MSAOBL_RSource.R")
Rcpp::sourceCpp(file="MSAOBL_CppSource.cpp", embeddedR = FALSE)
MSAOBL <- function(optFunction, popSize = 25, maxGeneration = 100, stopValue, maxWalkStep = 0.5, OBLStopPercent = 30, upperBounds, lowerBounds, numThreads = 1)
{
  if(length(lowerBounds) != length(upperBounds))
  {
    print("Bounds vectors of not same size");
    return (NULL)
  }
  
  for(i in 1:length(lowerBounds))
  {
    if(is.na(lowerBounds[i]))
    {
      lowerBounds[i] = -double.xmax
    }
    
    if(is.na(upperBounds[i]))
    {
      upperBounds[i] = double.xmax
    }
    if(lowerBounds[i] > upperBounds[i])
    {
      print("Lower Bound at index " + i + " higher than upper bound")
    }
  }
  
  
  
  if(is_xpt(optFunction) || numThreads == 1)
  {
    MothSearchOptCpp(optFunction, popSize, maxGeneration, stopValue, maxWalkStep, OBLStopPercent, upperBounds, lowerBounds, numThreads)
  }
  else
  {
    return (MothSearchOpt(optFunction, popSize, maxGeneration, stopValue, maxWalkStep, OBLStopPercent, upperBounds, lowerBounds, numThreads))
  }
}

optFunc <- function(x)
{
  return (rosenbrock(x))
}
optFuncCpp <- cppXPtr("double func(std::vector<double> x){
return (-20 * exp(-0.2 * sqrt(0.5 * (x[0] * x[0] + x[1] * x[1]))) -
		exp(0.5 * (cos(2 * M_PI*x[0]) + cos(2 * M_PI*x[1]))) + M_E + 20);
}")
result <- MSAOBL(optFunc, 25, 50, 0.0, 0.5, 0, rep(5, 2), rep(-5.0, 2), 1)
x <- 1:51
y <- result["BM_Fitness"]
y_vec <- y[[1]]
plot(x,y_vec,type="l")

