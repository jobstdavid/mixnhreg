// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
#include <boost/math/tools/roots.hpp>

//////////////////////////////
// CONSTANTS
//////////////////////////////
constexpr double PI = 3.14159265358979323846;
const double inv_sqrt_PI = 1/std::sqrt(PI);

//////////////////////////////
// HELPER FUNCTIONS BISECTION
//////////////////////////////

double pdist_root(const double& x,
                  const double& p,
                  const Eigen::VectorXd& location,
                  const Eigen::VectorXd& scale,
                  const Eigen::VectorXd& weight) {

  int K = location.size();
  double out = -p;

  for (int j = 0; j < K; ++j) {
    out += weight(j)*R::pnorm5(x, location(j), scale(j), 1, 0);
  }

  return(out);

}

// [[Rcpp::export]]
double root_finder(const Eigen::VectorXd& location,
                   const Eigen::VectorXd& scale,
                   const Eigen::VectorXd& weight,
                   const double& p,
                   double& xmin,
                   double& xmax,
                   const double& tol) {

  double out;

  if ((p < 0) | (p > 1)) {

    out = R_NaN;

  } else {

    int iter, maxit;
    double f_min, f_max, ox, of;
    bool stopCriterion;
    Eigen::Array<bool, Eigen::Dynamic, 1> check_finite(2);
    Eigen::VectorXd delta(2);
    Eigen::VectorXd bounds(2);

    iter = 0;
    maxit = 100;
    f_min = pdist_root(xmin, p, location, scale, weight);
    f_max = pdist_root(xmax, p, location, scale, weight);
    stopCriterion = f_min*f_max > 0;
    delta(0) = 0.1*std::fmax(1e-4, std::abs(xmin));
    delta(1) = 0.1*std::fmax(1e-4, std::abs(xmax));
    bounds(0) = xmin;
    bounds(1) = xmax;

    while (stopCriterion & bounds.array().isFinite().any()) {

      check_finite = bounds.array().isFinite();

      if (check_finite(0)) {
        ox = xmin;
        of = f_min;
        xmin = xmin - delta(0);
        f_min = pdist_root(xmin, p, location, scale, weight);
      }

      if (check_finite(1)) {
        ox = xmax;
        of = f_max;
        xmax = xmax + delta(1);
        f_max = pdist_root(xmax, p, location, scale, weight);
      }

      delta = 2*delta;
      bounds(0) = xmin;
      bounds(1) = xmax;
      stopCriterion = f_min*f_max > 0;
      iter += 1;
      if (iter >= maxit) {
        break;
      }
    }

    if (stopCriterion) {
      out = NA_REAL;
    } else {
      auto x = boost::math::tools::bisect(
        [p, location, scale, weight](double x){ return pdist_root(x, p, location, scale, weight); },
        bounds(0),
        bounds(1),
        [tol](double a, double b){return abs(a-b) < tol;}
      );
      out = (x.first+x.second)/2;
    }

  }

  return(out);

}


//////////////////////////////
// MIXTURE NORMAL DISTRIBUTION
//////////////////////////////

// [[Rcpp::export]]
Eigen::VectorXd dmixnorm_cpp(const Eigen::VectorXd& x,
                             const Eigen::MatrixXd& location,
                             const Eigen::MatrixXd& scale,
                             const Eigen::MatrixXd& weight,
                             const double& log_p){

  int N = x.size();
  int K = location.cols();
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      out(i) += weight(i,j)*R::dnorm4(x(i), location(i,j), scale(i,j), 0);
    }
  }

  if (log_p == 1) {
    out = out.array().log();
  }

  return(out);
}

// [[Rcpp::export]]
Eigen::VectorXd pmixnorm_cpp(const Eigen::VectorXd& q,
                             const Eigen::MatrixXd& location,
                             const Eigen::MatrixXd& scale,
                             const Eigen::MatrixXd& weight,
                             const double& log_p,
                             const double& lower_tail){

  int N = q.size();
  int K = location.cols();
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      out(i) += weight(i, j)*R::pnorm5(q(i), location(i, j), scale(i, j), 1, 0);
    }
  }

  if (lower_tail == 0) {
    out = 1-out.array();
  }

  if (log_p == 1) {
    out = out.array().log();
  }

  return(out);
}


// [[Rcpp::export]]
Eigen::VectorXd qmixnorm_cpp(const Eigen::VectorXd& p,
                             const Eigen::MatrixXd& location,
                             const Eigen::MatrixXd& scale,
                             const Eigen::MatrixXd& weight,
                             const double& log_p,
                             const double& lower_tail){

  int N = p.size();
  double xmin_norm, xmax_norm;
  double tol_norm = 1e-10;
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);

  if (log_p == 0) {
    if (lower_tail == 1) {
      for (int i = 0; i < N; ++i) {
        if (p(i) == 0 || p(i) == 1) {
          out(i) = pow(-1, 1+p(i))*R_PosInf;
        } else {
          xmin_norm = (location.row(i)-3*scale.row(i)).array().minCoeff();
          xmax_norm = (location.row(i)+3*scale.row(i)).array().maxCoeff();
          out(i) = root_finder(location.row(i), scale.row(i), weight.row(i), p(i),
              xmin_norm, xmax_norm, tol_norm);
        }
      }
    } else {
      for (int i = 0; i < N; ++i) {
        if (p(i) == 0 || p(i) == 1) {
          out(i) = pow(-1, p(i))*R_PosInf;
        } else {
          xmin_norm = (location.row(i)-3*scale.row(i)).array().minCoeff();
          xmax_norm = (location.row(i)+3*scale.row(i)).array().maxCoeff();
          out(i) = root_finder(location.row(i), scale.row(i), weight.row(i), 1-p(i),
              xmin_norm, xmax_norm, tol_norm);
        }
      }
    }
  } else {
    if (lower_tail == 1) {
      for (int i = 0; i < N; ++i) {
        if (p(i) == 0 || p(i) == 1) {
          out(i) = pow(-1, 1+p(i))*R_PosInf;
        } else {
          xmin_norm = (location.row(i)-3*scale.row(i)).array().minCoeff();
          xmax_norm = (location.row(i)+3*scale.row(i)).array().maxCoeff();
          out(i) = root_finder(location.row(i), scale.row(i), weight.row(i), exp(p(i)),
              xmin_norm, xmax_norm, tol_norm);
        }
      }
    } else {
      for (int i = 0; i < N; ++i) {
        if (p(i) == 0 || p(i) == 1) {
          out(i) = pow(-1, p(i))*R_PosInf;
        } else {
          xmin_norm = (location.row(i)-3*scale.row(i)).array().minCoeff();
          xmax_norm = (location.row(i)+3*scale.row(i)).array().maxCoeff();
          out(i) = root_finder(location.row(i), scale.row(i), weight.row(i), exp(1-p(i)),
              xmin_norm, xmax_norm, tol_norm);
        }
      }
    }
  }

  return(out);

}

// [[Rcpp::export]]
Eigen::VectorXd rmixnorm_cpp(const Eigen::VectorXd& p,
                             const Eigen::MatrixXd& location,
                             const Eigen::MatrixXd& scale,
                             const Eigen::MatrixXd& weight){

  int N = p.size();
  int K = location.cols();
  double cs;
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);

  for (int i = 0; i < N; ++i) {
    cs = 0;
    for (int j = 0; j < K; ++j) {
      cs += weight(i,j);
      if (p(i) < cs) {
        out(i) = R::qnorm5(p(i), location(i, j), scale(i, j), 1, 0);
        break;
      }
    }
  }

  return(out);

}



//////////////////////////////////////
// SCORES MIXTURE NORMAL DISTRIBUTION
/////////////////////////////////////

// [[Rcpp::export]]
Eigen::VectorXd logs_mixnorm_cpp(const Eigen::VectorXd& x,
                                 const Eigen::MatrixXd& location,
                                 const Eigen::MatrixXd& scale,
                                 const Eigen::MatrixXd& weight){

  Eigen::VectorXd out = -dmixnorm_cpp(x, location, scale, weight, 1);

  return(out);
}


double crps_mixnorm_helper_cpp(const double& location,
                               const double& scale){

  double out = location*(2*R::pnorm5(location/scale, 0, 1, 1, 0) - 1) + 2*scale*R::dnorm4(location/scale, 0, 1, 0);

  return(out);
}

// [[Rcpp::export]]
Eigen::VectorXd crps_mixnorm_cpp(const Eigen::VectorXd& x,
                                 const Eigen::MatrixXd& location,
                                 const Eigen::MatrixXd& scale,
                                 const Eigen::MatrixXd& weight){

  int N = x.size();
  int K = location.cols();
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);

  for (int i = 0; i < N; ++i) {

    for (int j = 0; j < K; ++j) {

      out(i) += weight(i, j)*crps_mixnorm_helper_cpp(x(i)-location(i, j), scale(i, j));
      out(i) -= inv_sqrt_PI*scale(i, j)*std::pow(weight(i, j), 2);

      for (int k = 0; k < j; ++k) {
        out(i) -= weight(i, j)*weight(i, k)*crps_mixnorm_helper_cpp(location(i, j)-location(i, k), std::sqrt(pow(scale(i, j), 2) + pow(scale(i, k), 2)));
      }

    }

  }

  return(out);

}


///////////////////////////////////////////
// GRADIENTS MIXTURE NORMAL DISTRIBUTION
//////////////////////////////////////////


// [[Rcpp::export]]
List grad_logs_mixnorm_cpp(const Eigen::VectorXd& x,
                           const Eigen::MatrixXd& location,
                           const Eigen::MatrixXd& scale,
                           const Eigen::MatrixXd& weight){

  int N = x.size();
  int K = location.cols();
  Eigen::VectorXd d = Eigen::VectorXd::Zero(K);
  Eigen::VectorXd pi = Eigen::VectorXd::Zero(K);
  Eigen::MatrixXd grad_location = Eigen::MatrixXd::Zero(N, K);
  Eigen::MatrixXd grad_scale = Eigen::MatrixXd::Zero(N, K);
  Eigen::MatrixXd grad_weight = Eigen::MatrixXd::Zero(N, K);

  for (int i = 0; i < N; ++i) {

    for (int j = 0; j < K; ++j) {
      d(j) = weight(i, j)*R::dnorm4(x(i), location(i, j), scale(i, j), 0);
    }
    pi = d/d.sum();

    for (int j = 0; j < K; ++j) {
      grad_location(i, j) = pi(j)*(location(i, j)-x(i))/pow(scale(i, j), 2);
      grad_scale(i, j) = pi(j)*(1 - pow((x(i)-location(i, j))/scale(i, j), 2));
      grad_weight(i, j) = weight(i, j) - pi(j);
    }

  }

  return List::create(Named("location") = grad_location, _["scale"] = grad_scale, _["weight"] = grad_weight);

}


// [[Rcpp::export]]
List grad_crps_mixnorm_cpp(const Eigen::VectorXd& x,
                           const Eigen::MatrixXd& location,
                           const Eigen::MatrixXd& scale,
                           const Eigen::MatrixXd& weight){

  int N = x.size();
  int K = location.cols();
  double CRPS, tmp;
  Eigen::MatrixXd grad_location = Eigen::MatrixXd::Zero(N, K);
  Eigen::MatrixXd grad_scale = Eigen::MatrixXd::Zero(N, K);
  Eigen::MatrixXd grad_weight = Eigen::MatrixXd::Zero(N, K);

  for (int i = 0; i < N; ++i) {

    CRPS = 0;
    for (int j = 0; j < K; ++j) {

      grad_location(i, j) = weight(i, j)*(1 - 2*R::pnorm5((x(i) - location(i, j))/scale(i, j), 0, 1, 1, 0));
      grad_scale(i, j) = 2*weight(i, j)*R::dnorm4((x(i)-location(i, j))/scale(i, j), 0, 1, 0) - inv_sqrt_PI*pow(weight(i, j), 2);
      grad_weight(i, j) = weight(i, j)*crps_mixnorm_helper_cpp(x(i)-location(i, j), scale(i, j)) - 2*inv_sqrt_PI*scale(i, j)*std::pow(weight(i, j), 2);
      CRPS -= grad_weight(i, j);

      for (int k = 0; k < j; k++) {
        grad_location(i, j) -= weight(i, k)*weight(i, j)*(1 - 2*R::pnorm5((location(i, k)-location(i, j))/std::sqrt(pow(scale(i, k), 2) + pow(scale(i, j), 2)), 0, 1, 1, 0));
        grad_scale(i, j) -= 2*weight(i, k)*weight(i, j)*scale(i, j)/std::sqrt(pow(scale(i, k), 2) + pow(scale(i, j), 2))*R::dnorm4((location(i, k)-location(i, j))/std::sqrt(pow(scale(i, k), 2) + pow(scale(i, j), 2)), 0, 1, 0);
        tmp = weight(i, k)*weight(i, j)*crps_mixnorm_helper_cpp(location(i, k)-location(i, j), std::sqrt(pow(scale(i, k), 2) + pow(scale(i, j), 2)));
        grad_weight(i, j) -= tmp;
        CRPS += 2*tmp;
      }

      for (int k = j + 1; k < K; k++) {
        grad_location(i, j) += weight(i, j)*weight(i, k)*(1 - 2*R::pnorm5((location(i, j)-location(i, k))/std::sqrt(pow(scale(i, j), 2) + pow(scale(i, k), 2)), 0, 1, 1, 0));
        grad_scale(i, j) -= 2*weight(i, j)*weight(i, k)*scale(i, j)/std::sqrt(pow(scale(i, j), 2) + pow(scale(i, k), 2))*R::dnorm4((location(i, j)-location(i, k))/std::sqrt(pow(scale(i, j), 2) + pow(scale(i, k), 2)), 0, 1, 0);
        grad_weight(i, j) -= weight(i, j)*weight(i, k)*crps_mixnorm_helper_cpp(location(i, j)-location(i, k), std::sqrt(pow(scale(i, j), 2) + pow(scale(i, k), 2)));
      }

      grad_scale(i, j) *= scale(i, j);

    }

    grad_weight.row(i) += weight.row(i)*CRPS;

  }

  return List::create(Named("location") = grad_location, _["scale"] = grad_scale, _["weight"] = grad_weight);

}

