library(foreach)
library(doParallel)

levyFlightDist <- function(beta, s)
{
  return (((beta - 1)*gamma(beta - 1)*sin((pi*(beta - 1)) / 2)) / pi * s^beta)
}

levyFlightMove <- function(maxWalkStep, currGen, coords, upperBounds, lowerBounds)
{
  levy = levyFlightDist(1.5, 2.0) # s = 2, zistit dobru hodnotu
  alpha = maxWalkStep / currGen^2
  
  resultCoords <- coords + levy * alpha
  
  for(i in 1:(length(resultCoords)))
  {
    if(resultCoords[i] < lowerBounds[i])
    {
      resultCoords[i] <- lowerBounds[i]
    }
    if(resultCoords[i] > upperBounds[i])
    {
      resultCoords[i] <- upperBounds[i]
    }
  }
  
  return (resultCoords)
}

EqFiveSixMove <- function(scaleFactor, goldenRatio, bestMothCoords, coords, upperBounds, lowerBounds)
{
  resultCoords <- scaleFactor * coords + goldenRatio * (bestMothCoords - coords)
  
  for(i in 1:(length(resultCoords)))
  {
    if(resultCoords[i] < lowerBounds[i])
    {
      resultCoords[i] <- lowerBounds[i]
    }
    if(resultCoords[i] > upperBounds[i])
    {
      resultCoords[i] <- upperBounds[i]
    }
  }
  
  return (resultCoords);
}

bestMothOBL <- function(bestMoth, coords, upperBounds, lowerBounds, optFunc)
{
  newMothCoords <- 2 * bestMoth - coords
  for(i in 1:(length(newMothCoords) - 1))
  {
    
    if ((newMothCoords[i] < lowerBounds[i]) || (newMothCoords[i] > upperBounds[i]))
    {
      newMothCoords[i] = runif(1, lowerBounds[i],  upperBounds[i])
    }
  }
  newMothCoords[length(newMothCoords)] <- optFunc(newMothCoords[1:(length(newMothCoords)- 1)]);
  
  return (newMothCoords);
}

boundsOBL <- function(upperBounds, lowerBounds, optFunc, coords)
{
  #upperbounds and lowerBounds need added zero (any number) at the end to match length of moth vector
  newMothCoords <- upperBounds + lowerBounds - coords
  
  newMothCoords[length(newMothCoords)] <- optFunc(newMothCoords[1:(length(newMothCoords) - 1)])
  
  return (newMothCoords)
  
}

MothSearchOpt <- function(optFunction, popSize, maxGeneration, stopValue, maxWalkStep, OBLStopPercent, upperBounds, lowerBounds, numThreads)
{
  #parallel inits
  exports <- c("levyFlightMove", "levyFlightDist", "boundsOBL", "bestMothOBL", "EqFiveSixMove")
  cl <- parallel::makeCluster(numThreads)
  doParallel::registerDoParallel(cl)
  
  #constant inits
  goldenRatio <- (1+sqrt(5))/2
  scaleFactor <- rnorm(1, 0.5, 0.1)
  OBLStop <- maxGeneration * OBLStopPercent / 100
  
  #population init
  popMatrix <- matrix(data = NA, nrow = length(lowerBounds) + 1, ncol = popSize)
  for(i in 1:length(lowerBounds))
  {
    
    popMatrix[i,] <- runif(n=popSize, min=lowerBounds[i], max=upperBounds[i])
  }
  popMatrix[length(lowerBounds) + 1,] <- numeric(popSize)
  popMatrix <- t(popMatrix)
  
  #compute fitness of initial pop
  for(i in 1:nrow(popMatrix))
  {
    popMatrix[i,ncol(popMatrix)] <- optFunction(popMatrix[i,1:ncol(popMatrix) - 1])
  }
  #added zero to bounds to match length of a Moth representations
  lowerBounds <- c(lowerBounds, 0.0)
  upperBounds <- c(upperBounds, 0.0)
  dimensions <- length(lowerBounds)
  
  #output inits
  populationFitnesses <- matrix(data=NA, nrow = maxGeneration + 1, ncol=popSize)
  populationFitnesses[1,] <- popMatrix[,dimensions]
  populationCoords <- vector(mode = "list", length = maxGeneration + 1)
  populationCoords[[1]] <- popMatrix[,1:(dimensions - 1)]
  bestMothFitnesses <- numeric(maxGeneration + 1)
  bestMothFitnesses[1] <- popMatrix[1,dimensions]
  bestMothCoords <- matrix(data=NA, nrow = maxGeneration + 1, ncol=dimensions - 1)
  bestMothCoords[1,] <- popMatrix[1,1:(dimensions - 1)]
  
  for(j in 1:maxGeneration)
  {
    currGen <- j
    #sort pop by fitness
    popMatrix <- popMatrix[order(popMatrix[, ncol(popMatrix)]),]
    bestMoth <- popMatrix[1,]
    popMatrix <- foreach(i=1:popSize, .combine=rbind, .export=exports) %dopar%
    {
      # better half of population move (exploration)
      if(i <= popSize/2)
      {
        newPositionMoth <- levyFlightMove(maxWalkStep, currGen, popMatrix[i,], upperBounds, lowerBounds)
        newPositionMoth[dimensions] <- optFunction(newPositionMoth[1:(dimensions - 1)])
        if(currGen < OBLStop){
          OBLMoth <- bestMothOBL(bestMoth, popMatrix[i,], upperBounds, lowerBounds, optFunction)
          if (OBLMoth[dimensions] < newPositionMoth[dimensions])
          {
            return (OBLMoth)
          }
          else
          {
            return (newPositionMoth)
          }
        }
        else
        {
          return (newPositionMoth)
        }
      }
      else #worse half of population move (exploitation)
      {
        randDecider <- runif(1, 0.0, 1.0);
        
        if (randDecider < 0.5)
        {
          newPositionMoth <- EqFiveSixMove(scaleFactor, goldenRatio, bestMoth, popMatrix[i,], upperBounds, lowerBounds)
        }
        else
        {
          newPositionMoth <- EqFiveSixMove(scaleFactor, 1 / goldenRatio, bestMoth, popMatrix[i,], upperBounds, lowerBounds)
        }
        newPositionMoth[dimensions] <- optFunction(newPositionMoth[1:(dimensions - 1)])
        if(currGen < OBLStop)
        {
          OBLMoth <- boundsOBL(upperBounds, lowerBounds, optFunction, popMatrix[i,])
          if (OBLMoth[dimensions] < newPositionMoth[dimensions])
          {
            return (OBLMoth)
          }
          else
          {
            return (newPositionMoth)
          }
        }
        else return (newPositionMoth)
      }
    }
    
    if(bestMoth[dimensions] <= stopValue)
    {
      break
    }
    
    #writing output values
    populationFitnesses[currGen + 1,] <- popMatrix[,dimensions]
    populationCoords[[currGen + 1]] <- popMatrix[,1:(dimensions - 1)]
    bestMothFitnesses[currGen + 1] <- bestMoth[dimensions]
    bestMothCoords[currGen + 1,] <- bestMoth[1:(dimensions - 1)]
    
  }
  parallel::stopCluster(cl)
  
  resultList <- list(bestMothFitnesses, bestMothCoords, populationFitnesses, populationCoords)
  names(resultList) <- c("BM_Fitness","BM_Coords", "POP_Fitness", "POP_Coords")
  return(resultList)
}