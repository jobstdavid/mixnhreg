// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

//////////////////////////////
// CONSTANTS
//////////////////////////////
constexpr double PI = 3.14159265358979323846;

///////////////////////////////////////////
// HELPER FUNCTIONS OPTIM
//////////////////////////////////////////

// [[Rcpp::export]]
Eigen::VectorXd unlist_cpp(const Rcpp::List& list) {
  std::size_t n = list.size();

  // Calculate the total length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i) {
    total_length += Rf_length(list[i]);
  }

  // Prepare an Eigen::VectorXd for the output
  Eigen::VectorXd output(total_length);

  // Loop through the list and fill the Eigen vector
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i) {
    // Extract each element as an Rcpp::NumericVector
    Rcpp::NumericVector el = list(i);

    // Map the Rcpp vector to an Eigen::Map for efficient data access
    Eigen::Map<Eigen::VectorXd> eigen_el(el.begin(), el.size());

    // Place this segment into the output Eigen vector
    output.segment(index, el.size()) = eigen_el;

    // Update the index
    index += el.size();
  }

  // Return the Eigen vector
  return output;
}

// [[Rcpp::export]]
Eigen::VectorXd grad_optim_cpp(const Rcpp::List& X,
                               const Eigen::MatrixXd& G){

  // X: list of predictor matrices
  // G: matrix of gradients

  int D = G.cols();
  Rcpp::List output;

  for (int i = 0; i < D; ++i) {
    Eigen::MatrixXd predictor(Rcpp::as<Eigen::MatrixXd>(X.at(i)));
    output.push_back(predictor.transpose() * G.col(i));
  }

  return(unlist_cpp(output));

}

// [[Rcpp::export]]
Eigen::MatrixXd get_parameter_cpp(const Eigen::MatrixXd& X,
                                  const std::string& link,
                                  const bool& inverse) {

  // X: predictor matrix
  // link: (inverse) link function

  Eigen::MatrixXd output;

  // get parameter based on (inverse) link function
  if (link == "id") {
    output = X;
  } else if (link == "log") {
    if (!inverse) {
      output = X.array().log();
    } else {
      output = X.array().exp();
    }
  } else if (link == "sqrt") {
    if (!inverse) {
      output = X.array().sqrt();
    } else {
      output = X.array().pow(2);
    }
  } else if (link == "logit") {
    if (!inverse) {
      output = (X.array()/(1-X.array())).log();
    } else {
      output = 1/(1+(-X.array()).exp());
    }
  } else if (link == "loglog") {
    if (!inverse) {
      output = (-X.array().log()).log();
    } else {
      output = (-X.array().exp()).exp();
    }
  } else if (link == "cloglog") {
    if (!inverse) {
      output = (-(1-X.array()).log()).log();
    } else {
      output = 1-((-X.array().exp()).exp());
    }
  } else if (link == "cauchit") {
    if (!inverse) {
      output = (PI*(X.array() - 1/2)).tan();
    } else {
      output = X.array().atan()/PI + 1/2;
    }
  } else if (link == "atanh") {
    if (!inverse) {
      output = X.array().atanh();
    } else {
      output = X.array().tanh();
    }
  } else if (link == "softmax") {
    output = X.array().exp();
    Eigen::VectorXd s = output.rowwise().sum();
    output = output.array().colwise()/s.array();
  }

  return(output);

}

// [[Rcpp::export]]
Rcpp::List get_parameter_optim_cpp(const Rcpp::List& X,
                                   const Eigen::VectorXd& x,
                                   const std::vector<std::string>& links,
                                   const int& N,
                                   const int& K){

  // X: list of predictor matrices
  // x: vector of coefficients

  int l = 0, iter = 0;
  int u = -1;
  Rcpp::List output;

  for (int j = 0; j < 3; ++j) {

    // get predictor
    Eigen::MatrixXd predictor(N, K);
    for (int i = 0; i < K; ++i) {
      Eigen::MatrixXd Y(Rcpp::as<Eigen::MatrixXd>(X.at(iter)));
      u += Y.cols();
      predictor.col(i) = Y * x.segment(l, u-l+1);
      l = u+1;
      iter += 1;
    }

    // get parameter based on link function
    output.push_back(get_parameter_cpp(predictor, links[j], TRUE));

  }

  return(output);

}

///////////////////////////////////////////
// HELPER FUNCTIONS GRADIENT-BOOSTING
//////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List get_parameter_boost_cpp(const Rcpp::List& X,
                                   const Rcpp::List& x,
                                   const std::vector<std::string>& links,
                                   const int& N,
                                   const int& K){

  // X: list of predictor matrices
  // x: list of coefficients of where length(X) = length(x)

  int iter = 0;
  Rcpp::List output;

  for (int j = 0; j < 3; ++j) {

    // get predictor
    Eigen::MatrixXd predictor(N, K);
    for (int i = 0; i < K; ++i) {
      Eigen::MatrixXd Y(Rcpp::as<Eigen::MatrixXd>(X.at(iter)));
      Eigen::VectorXd p(Rcpp::as<Eigen::VectorXd>(x.at(iter)));
      predictor.col(i) = Y * p;
      iter += 1;
    }

    // get parameter based on link function
    output.push_back(get_parameter_cpp(predictor, links[j], TRUE));

  }

  return(output);

}

// [[Rcpp::export]]
Rcpp::List get_basefits_cpp(const Rcpp::List& X,
                            const Rcpp::List& G,
                            const int& N,
                            const int& K){

  // X: list of predictor matrices
  // G: list of gradient matrices

  Rcpp::List output;
  int iter = 0;

  for (int i = 0; i < 3; ++i) { // iterations of location, scale, weight parameter
    Eigen::MatrixXd Y(Rcpp::as<Eigen::MatrixXd>(G.at(i)));
    for (int j = 0; j < K; ++j) { // iterates over number of components
      Eigen::MatrixXd Z(Rcpp::as<Eigen::MatrixXd>(X.at(iter)));
      output.push_back((Z.transpose() * (-1*Y.col(j)))/N);
      iter += 1;
    }
  }

  return(output);

}

// [[Rcpp::export]]
Rcpp::List update_coef_cpp(const Eigen::VectorXd& y,
                           const Rcpp::List& X,
                           const Rcpp::List& x,
                           const Rcpp::List& B,
                           const std::vector<std::string>& links,
                           const Rcpp::Function loss_fun,
                           const double& nu,
                           const int& N,
                           const int& K){


  // B: list of basefits
  // x list of latest coefficients
  int index;
  int D = 3*K;
  Rcpp::List Z, tmp, output(D);
  Eigen::VectorXd bf, coef, loss;
  Eigen::VectorXd scores(D);

  for (int j = 0; j < D; ++j) {

    tmp = clone(x);
    coef = tmp(j);
    bf = B(j);
    bf.array().abs().maxCoeff(&index);
    coef(index) += nu*bf(index);
    tmp(j) = coef;
    output(j) = tmp;

    Z = get_parameter_boost_cpp(X, tmp, links, N, K);

    loss = Rcpp::as<Eigen::VectorXd>(loss_fun(y, Z[0], Z[1], Z[2]));
    scores(j) = loss.sum();

  }

  scores.minCoeff(&index);
  output = output(index);

  return(output);

}

// [[Rcpp::export]]
Rcpp::List boost_noncylcic_cpp(const Eigen::VectorXd& y,
                               const Rcpp::List& X,
                               Rcpp::List& x,
                               const std::vector<std::string>& links,
                               const Rcpp::Function loss_fun,
                               const Rcpp::Function grad_fun,
                               const Rcpp::Function nll_fun,
                               const double& nu,
                               const int& maxit,
                               const int& N,
                               const int& K){

  Rcpp::List Z;
  Eigen::VectorXd x0 = unlist_cpp(x);
  Eigen::VectorXd nll, nll_path(maxit+1);
  Eigen::MatrixXd coef_path(maxit+1, x0.size());
  coef_path.row(0) = x0;

  for (int i = 0; i < maxit; ++i) {

    Z = get_parameter_boost_cpp(X, x, links, N, K);
    nll = Rcpp::as<Eigen::VectorXd>(nll_fun(y, Z[0], Z[1], Z[2]));
    nll_path(i) = nll.sum();

    Z = grad_fun(y, Z[0], Z[1], Z[2]);

    Z = get_basefits_cpp(X, Z, N, K);

    x = update_coef_cpp(y, X, x, Z, links, loss_fun, nu, N, K);

    coef_path.row(i+1) = unlist_cpp(x);

  }

  Z = get_parameter_boost_cpp(X, x, links, N, K);
  nll = Rcpp::as<Eigen::VectorXd>(nll_fun(y, Z[0], Z[1], Z[2]));
  nll_path(maxit) = nll.sum();

  Rcpp::List output = List::create(_["coef_path"] = coef_path,
                                   _["nll_path"] = nll_path);

  return(output);

}

// [[Rcpp::export]]
Eigen::MatrixXd get_loss_cv_cpp(const Eigen::VectorXd& y,
                                const Rcpp::List& X,
                                const Eigen::MatrixXd& x,
                                const std::vector<std::string>& links,
                                const Rcpp::Function loss_fun,
                                const int& maxit,
                                const int& N,
                                const int& K){

  // Rcpp::List output;
  Eigen::MatrixXd output(N, maxit+1);
  Eigen::VectorXd loss;
  Rcpp::List parameter;

  for (int m = 0; m < maxit+1; ++m) {

    parameter = get_parameter_optim_cpp(X, x.row(m), links, N, K);
    loss = Rcpp::as<Eigen::VectorXd>(loss_fun(y, parameter[0], parameter[1], parameter[2]));
    output.col(m) = loss;

  }

  return(output);

}





