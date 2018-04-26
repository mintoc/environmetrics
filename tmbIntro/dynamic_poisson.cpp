#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  // DATA
  DATA_VECTOR(y);
  DATA_VECTOR(ypres);
  PARAMETER(lnsigma);
  PARAMETER_VECTOR(lnlambda);
  // PRELIMINARY CALCULATIONS
  int n = y.size(); // number of observations
  Type sigma = exp(lnsigma);
  // PROCEDURE
  Type nll = 0.0; // initialize negative log likelihood
  // likelihood of the random effects
  for(int i = 1; i < n; i++){ 
    nll -= dnorm(lnlambda(i), lnlambda(i-1), sigma, true);
  }
  vector<Type> lambda = exp(lnlambda);
  // loop over the observations
  for(int i = 0; i < n; i++){
    if(ypres(i) > 0){
      nll -= y(i) * lnlambda(i) - lambda(i) -lfactorial(y(i));
    }
  }
  return nll;
}
