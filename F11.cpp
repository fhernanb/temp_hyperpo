#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <complex>

// Global variable to cache the result
double f11_cache;

// Function to check the stopping condition
bool stopping(double diff, double tol) {
  return std::abs(diff) < tol;
}

// [[Rcpp::export]]
double f11_cpp(double gamma, double lambda, int maxiter_series = 10000, double tol = 1.0e-10) {
  double fac  = 1.0;
  double temp = 1.0;
  double L    = gamma;
  double series = temp;
  
  for (int n = 1; n <= maxiter_series; ++n) {
    fac = fac * lambda / L;
    series = temp + fac;
    
    if (stopping(series - temp, tol)) {
      f11_cache = series;  // Assuming series is already real in this context
      return f11_cache;
    }
    
    temp = series;
    L += 1;
  }
  
  if (tol >= 0) {
    std::cerr << "Warning: Tolerance is not met" << std::endl;
  }
  
  f11_cache = series;  // Assuming series is already real in this context
  return f11_cache;
}


// Vectorized version for dHYPERPO
// [[Rcpp::export]]
double dHYPERPO_cpp(double x, double mu=1, double sigma=1, bool log=false) {
  if (sigma <= 0 || mu <= 0) {
    throw std::runtime_error("parameter sigma and mu must be positive!");
  }
  
  double res;
  
  if (x < 0) {
    res = std::log(0);
  }
  else {
    double p1 = x * std::log(mu) - lgamma(sigma + x) + lgamma(sigma);
    double temp_f11 = f11_cpp(sigma, mu);
    double p2 = std::log(temp_f11);
    res = p1 - p2;
  }
  
  if (log) {
    return res;
  } else {
    return std::exp(res);
  }
}

// [[Rcpp::export]]
NumericVector dHYPERPO_vec(NumericVector x, NumericVector mu, 
                           NumericVector sigma, LogicalVector log) {
  int n = x.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; ++i) {
    out[i] = dHYPERPO_cpp(x[i], mu[i], sigma[i], log[i]);
  }
  
  return out;
}


// Function to be used inside obtaining_lambda
// [[Rcpp::export]]
double fun(double x, double gamma, double media) {
  return x - (gamma - 1) * (f11_cpp(gamma, x) - 1) / f11_cpp(gamma, x) - media;
}

// Bisection method to find root, similar to uniroot in R
double find_root(std::function<double(double)> f, double lower, double upper, double tol = 1.0e-10, int max_iter = 1000) {
  double mid;
  for (int iter = 0; iter < max_iter; ++iter) {
    mid = (lower + upper) / 2.0;
    double f_mid = f(mid);
    
    if (std::abs(f_mid) < tol) {
      return mid;
    }
    
    if (f(lower) * f_mid < 0) {
      upper = mid;
    } else {
      lower = mid;
    }
  }
  
  //throw std::runtime_error("Tolerance not met in root finding");
  return 1;
}

// [[Rcpp::export]]
double obtaining_lambda_cpp(double media, double gamma) {
  double result;
  auto fun_wrapper = [&](double x) { return fun(x, gamma, media); };
  
  if (gamma == 1) {
    result = media;
  } else {
    double lower = std::min(media, std::max(media + gamma - 1, gamma * media));
    double upper = std::max(media, std::min(media + gamma - 1, gamma * media));
    result = find_root(fun_wrapper, lower, upper);
  }
  
  return result;
}


