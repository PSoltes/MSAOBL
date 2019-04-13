// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppParallel)]]

#define _USE_MATH_DEFINES
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>
#include <Rcpp.h>
#include <omp.h>
#include "functionWrappers.h"

double levyFlightDist(double beta, double s)
{
  return ((beta - 1)*tgamma(beta - 1)*sin((M_PI*(beta - 1)) / 2)) / M_PI * pow(s, beta);
}

double randomFromInterval(double max, double min)
{
  return (max - min) * ((double)rand() / (double)RAND_MAX) + min;
}

class Moth
{
public:
  Moth() = default;
  
  Moth(std::vector<double> newMothCoords)
    :_coords(newMothCoords)
  {
    _fitness = 0;
  }
  
  std::vector<double> getCoords()
  {
    return this->_coords;
  }
  
  double getFitness()
  {
    return this->_fitness;
  }
  
  void computeFitness(const std::unique_ptr<FunctionBase>& optFunction)
  {
    this->_fitness = optFunction->use(this->_coords);
  }
  
  void levyFlightMove(double maxWalkStep, size_t currGen, const std::vector<double>& upperBounds, const std::vector<double>& lowerBounds)
  {
    double levy = levyFlightDist(1.5, 2.0); // s = 2, zistit dobru hodnotu
    double alpha = maxWalkStep / pow(currGen, 2);
    
    for (size_t i = 0; i < this->_coords.size(); i++)
    {
      this->_coords[i] = round(this->_coords[i] + levy * alpha);
      
      if(this->_coords[i] < lowerBounds[i])
      {
        this->_coords[i] = lowerBounds[i];
      }
      if(this->_coords[i] > upperBounds[i])
      {
        this->_coords[i] = upperBounds[i];
      }
    }
    
    
  }
  
  void EqFiveSixMove(double scaleFactor, double goldenRatio, const Moth& bestMoth, const std::vector<double>& upperBounds, const std::vector<double>& lowerBounds)
  {
    for (size_t i = 0; i < this->_coords.size(); i++)
    {
      this->_coords[i] = round(scaleFactor * (this->_coords[i] + goldenRatio * (bestMoth._coords[i] - this->_coords[i])));
      
      if(this->_coords[i] < lowerBounds[i])
      {
        this->_coords[i] = lowerBounds[i];
      }
      if(this->_coords[i] > upperBounds[i])
      {
        this->_coords[i] = upperBounds[i];
      }
    }
    
  }
  
  const bool operator<(const Moth& rhs) const
  {
    return this->_fitness < rhs._fitness;
  }
  
  Moth bestMothOBL(const Moth& bestMoth, const std::vector<double>& upperBounds, const std::vector<double>& lowerBounds, const std::unique_ptr<FunctionBase>& optFunc)
  {
    std::vector<double> newMothCoords(this->_coords.size(), 0.0);
    
    for (size_t i = 0; i < this->_coords.size(); i++)
    {
      newMothCoords[i] = round(2 * bestMoth._coords[i] - this->_coords[i]);
      if(newMothCoords[i] < lowerBounds[i])
      {
        newMothCoords[i] = lowerBounds[i];
      }
      if(newMothCoords[i] > upperBounds[i])
      {
        newMothCoords[i] = upperBounds[i];
      }
    }
    
    Moth returnMoth(newMothCoords);
    
    returnMoth.computeFitness(optFunc);
    
    return returnMoth;
  }
  
  Moth boundsOBL(const std::vector<double>& upperBounds, const std::vector<double>& lowerBounds, const std::unique_ptr<FunctionBase>& optFunc)
  {
    std::vector<double> newMothCoords(this->_coords.size(), 0.0);
    
    for (size_t i = 0; i < this->_coords.size(); i++)
    {
      newMothCoords[i] = round(upperBounds[i] + lowerBounds[i] - this->_coords[i]);
    }
    
    Moth returnMoth(newMothCoords);
    
    returnMoth.computeFitness(optFunc);
    
    return returnMoth;
    
  }
  
private:
  std::vector<double> _coords;
  double _fitness;
};


std::vector<Moth> generatePop(size_t popSize, const std::vector<double>& upperBounds, const std::vector<double>& lowerBounds)
{
  std::vector<Moth> returnVector;
  size_t dimensions = upperBounds.size();
  for (size_t i = 0; i < popSize; i++)
  {
    Rcpp::NumericVector currentMothCoords(dimensions);
    for (size_t j = 0; j < dimensions; j++)
    {
      currentMothCoords[j] = R::runif(lowerBounds[j],upperBounds[j]);
    }
    
    returnVector.emplace_back(Rcpp::as<std::vector<double>>(currentMothCoords));
    
  }
  
  
  
  return returnVector;
}
// [[Rcpp::export]]
Rcpp::List MothSearchOptCppRounded(SEXP optFunc, size_t popSize, int maxGeneration, double stopValue, double maxWalkStep, int OBLStopPercent, Rcpp::NumericVector rUpperBounds, Rcpp::NumericVector rLowerBounds, int numThreads)
{
  size_t currentGeneration = 1;
  
  std::random_device rd{};
  std::mt19937 gen{ rd() };
  std::normal_distribution<double> distribution(0.5, 0.1);
  double goldenRatio = 1 + sqrt(5) / 2;
  double scaleFactor = distribution(gen);
  size_t OBLStop = maxGeneration * OBLStopPercent / 100;
  
  std::vector<double> upperBounds = Rcpp::as<std::vector<double>>(rUpperBounds);
  std::vector<double> lowerBounds = Rcpp::as<std::vector<double>>(rLowerBounds);
  Moth bestMoth;
  std::vector<Moth> population = generatePop(popSize, upperBounds, lowerBounds);
  std::vector<double> bestMothFitnesses;
  std::vector<std::vector<double>> bestMothCoords;
  std::vector<std::vector<double>> populationFitnesses;
  std::vector<std::vector<std::vector<double>>>  populationCoords;
  
  if(numThreads != 0)
  {
    omp_set_num_threads(numThreads);
  }
  
  
  bool parallel = false;
  std::unique_ptr<FunctionBase> _optFunc;
  
  if(TYPEOF(optFunc) == EXTPTRSXP)
  {
    parallel = true;
    Rcpp::XPtr<function> funcXPtr(optFunc);
    _optFunc = std::move(std::unique_ptr<FunctionCpp>(new FunctionCpp(funcXPtr)));
  }
  else
  {
    _optFunc = std::move(std::unique_ptr<FunctionR>(new FunctionR(optFunc)));
  }
  
  
  
#pragma omp parallel for if(parallel)
  for (int i = 0; i < population.size(); i++)
  {
    population[i].computeFitness(_optFunc);
  }
  //sort population by fitness
  std::sort(population.begin(), population.end());
  
  //output inits
  std::vector<double> thisGenFitnesess;
  std::vector<std::vector<double>> thisGenCoords;
  bestMothFitnesses.push_back(population[1].getFitness());
  bestMothCoords.push_back(population[1].getCoords());
  for(auto currMoth : population)
  {
    thisGenFitnesess.push_back(currMoth.getFitness());
    thisGenCoords.push_back(currMoth.getCoords());
  }
  populationFitnesses.push_back(thisGenFitnesess);
  populationCoords.push_back(thisGenCoords);
  
  for (; currentGeneration < maxGeneration; currentGeneration++)
  {
    std::sort(population.begin(), population.end());
    if(!(bestMoth < population[0]) || currentGeneration == 1)
    {
      bestMoth = population[0];
    }
#pragma omp parallel for if(parallel)
    for (int i = 0; i < population.size() / 2; i++)
    {
      population[i].levyFlightMove(maxWalkStep, currentGeneration, upperBounds, lowerBounds);
      population[i].computeFitness(_optFunc);
      if(currentGeneration < OBLStop)
      {
        try {
          Moth OBLMoth = population[i].bestMothOBL(bestMoth, upperBounds, lowerBounds, _optFunc);
          
          if (OBLMoth < population[i])
          {
            std::swap(OBLMoth, population[i]);
          }
        }
        catch (const std::exception& e)
        {
          std::terminate(); //to do something useful
        }
      }
      
    }
#pragma omp parallel for if(parallel)
    for (int j = population.size() / 2 + 1; j < population.size(); j++)
    {
      double randDecider = ((double)rand() / (RAND_MAX));
      
      if (randDecider < 0.5)
      {
        population[j].EqFiveSixMove(scaleFactor, goldenRatio, bestMoth, upperBounds, lowerBounds);
      }
      else
      {
        population[j].EqFiveSixMove(scaleFactor, 1 / goldenRatio, bestMoth, upperBounds, lowerBounds);
      }
      population[j].computeFitness(_optFunc);
      if(currentGeneration < OBLStop)
      {
        Moth OBLMoth = population[j].boundsOBL(upperBounds, lowerBounds, _optFunc);
        if (OBLMoth < population[j])
        {
          std::swap(OBLMoth, population[j]);
        }
      }
    }
    
    bestMothFitnesses.push_back(bestMoth.getFitness());
    bestMothCoords.push_back(bestMoth.getCoords());
    for(auto currMoth : population)
    {
      thisGenFitnesess.push_back(currMoth.getFitness());
      thisGenCoords.push_back(currMoth.getCoords());
    }
    populationFitnesses.push_back(thisGenFitnesess);
    populationCoords.push_back(thisGenCoords);
    
    if(bestMoth.getFitness() <= stopValue)
    {
      break;
    }
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("BM_Fitness") = bestMothFitnesses,
    Rcpp::Named("BM_Coords") = bestMothCoords,
    Rcpp::Named("POP_Fitness") = populationFitnesses,
    Rcpp::Named("POP_Coords")  = populationCoords
  );
}

