#ifndef R_FUNCTION_WRAPPERS
#define R_FUNCTION_WRAPPERS
#include <Rcpp.h>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

typedef double (*function)(std::vector<double>);


class FunctionBase
{
public:
  virtual double use(std::vector<double>) = 0;
private:
};

class FunctionCpp : public FunctionBase
{
public:
  FunctionCpp(Rcpp::XPtr<function> fn)
  {
    this->_fn = *fn;
  }
  
  double use(std::vector<double> x) override
  {
    return _fn(x);
  }
  
private:
  function _fn;
};

class FunctionR : public FunctionBase
{
public:
  FunctionR(SEXP fn) : _fn(fn)
  {}
  
  double use(std::vector<double> x) override
  {
    
    return Rcpp::as<double>(_fn(Rcpp::wrap(x)));
  }
  
private:
  Rcpp::Function _fn;
};
#endif
