#include <RcppEigen.h>
#include <iostream>
#include <thread>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace Eigen;

inline double relu(double g) {
  return (g > 1e-8 ? g : 0.0);
}

inline double sign(double x) {
  return (x > 1e-8) ? 1.0 : ((x < -1e-8) ? -1.0 : 0.0);
}

// Coordinate descent lasso with per-coordinate lambda (lambda_vec length P)
VectorXd lasso_weighted(const MatrixXd X, const VectorXd Y, const VectorXd B,
                        const VectorXd lambda_vec, int niter, int fixd = -1) {
  int P = X.cols();
  VectorXd E = Y - X * B;
  VectorXd B_next = B;
  for (int it = 0; it < niter; ++it) {
    for (int p = 0; p < P; ++p) {
      VectorXd x = X.col(p);
      if (fixd == p + 1) {
        B_next[p] = 1.0;
        E += x;
      } else {
        if (std::abs(B_next[p]) > 1e-8) {
          E += B_next[p] * x;
        }
        double mp = x.dot(E);
        double lam = (p < lambda_vec.size()) ? lambda_vec[p] : lambda_vec[lambda_vec.size() - 1];
        B_next[p] = sign(mp) * relu(std::abs(mp) - lam) / x.squaredNorm();
      }
      E -= B_next[p] * x;
    }
  }
  return B_next;
}

// [[Rcpp::export]]
VectorXd lasso(const MatrixXd X, const VectorXd Y, const VectorXd B,
               const double lambda, int niter, int fixd = -1) {
  int P = X.cols();
  VectorXd E = Y - X * B;
  VectorXd B_next = B;
  for (int i = 0; i < niter; i++) {
    for (int p = 0; p < P; p++) {
      VectorXd x = X.col(p);
      if (fixd == p + 1) {
        B_next[p] = 1.0;
        E += x;
      } else {
        if (std::abs(B_next[p]) > 1e-8) {
          E += B_next[p] * x;
        }
        double mp = x.dot(E);
        B_next[p] = sign(mp) * relu(std::abs(mp) - lambda) / x.squaredNorm();
      }
      E -= B_next[p] * x;
    }
  }
  return B_next;
}

// [[Rcpp::export]]
MatrixXd fit_U_VU_const(const Map<MatrixXd> WX, const Map<MatrixXd> W, const Map<MatrixXd> V,
                        const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double gamma = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1) {
  int D = U.cols();
  MatrixXd rVtV = rho * V.transpose() * V;
  MatrixXd rhs = WX + rho * V.transpose() - V.transpose() * theta;
  MatrixXd U_next = U;

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&, t]() {
      for (int d = t; d < D; d += nthreads) {
        MatrixXd A = (W.col(d).array() + gamma).matrix().asDiagonal();
        A += rVtV;
        U_next.col(d) = lasso(A, rhs.col(d), U.col(d), 0.0, solve_its, fixd ? d + 1 : -1);
      }
    });
  }
  for (auto& thread : threads) thread.join();
  return U_next;
}

// [[Rcpp::export]]
MatrixXd fit_U_UV_const(const Map<MatrixXd> WX, const Map<MatrixXd> W, const Map<MatrixXd> V,
                        const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double gamma = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1,
                        const bool adaptive_lambda = false, const NumericVector beta_obs = NumericVector()) {
  int D = U.cols();
  MatrixXd rVVt = rho * V * V.transpose();
  MatrixXd rhs = WX.transpose() + rho * V - V * theta.transpose();
  MatrixXd U_next = U;

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&, t]() {
      for (int d = t; d < D; d += nthreads) {
        MatrixXd A = (W.row(d).array() + gamma).matrix().asDiagonal();
        A += rVVt;
        U_next.row(d) = lasso(A, rhs.col(d), U.row(d).transpose(), 0.0, solve_its, fixd ? d + 1 : -1).transpose();
      }
    });
  }
  for (auto& thread : threads) thread.join();
  return U_next;
}

// [[Rcpp::export]]
MatrixXd fit_V_VU_const(const Map<MatrixXd> V, const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double lambda = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1,
                        const bool adaptive_lambda = false, const NumericVector beta_obs = NumericVector()) {
  int D = V.cols();
  MatrixXd V_next = V;
  MatrixXd glm_Y = MatrixXd::Identity(D, D) - (theta.transpose().array() / rho).matrix();

  std::vector<double> beta_vec(beta_obs.begin(), beta_obs.end());
  double beta_mean = 1.0;
  if (adaptive_lambda && beta_vec.size() == D) {
    double sum = 0.0;
    for (int i = 0; i < D; i++) sum += std::abs(beta_vec[i]);
    beta_mean = sum / D;
  }

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&, t]() {
      for (int d = t; d < D; d += nthreads) {
        double lambda_d = lambda;
        if (adaptive_lambda && beta_vec.size() == D) {
          lambda_d *= std::abs(beta_vec[d]) / (beta_mean + 1e-8);
        }
        // Use per-coordinate lambda (here uniform across coordinates for this row)
        VectorXd lambda_vec = VectorXd::Constant(D, lambda_d / rho);
        V_next.row(d) = lasso_weighted(U.transpose(), glm_Y.col(d), V.row(d).transpose(), lambda_vec, solve_its, fixd ? d + 1 : -1).transpose();
      }
    });
  }
  for (auto& thread : threads) thread.join();
  return V_next;
}

// [[Rcpp::export]]
MatrixXd fit_V_UV_const(const Map<MatrixXd> V, const Map<MatrixXd> U, const Map<MatrixXd> theta, const double rho,
                        const double lambda = 0, const int solve_its = 10,
                        const bool fixd = false, const int nthreads = 1,
                        const bool adaptive_lambda = false, const NumericVector beta_obs = NumericVector()) {
  int D = V.cols();
  MatrixXd V_next = V;
  MatrixXd glm_Y = MatrixXd::Identity(D, D) - (theta.array() / rho).matrix();

  std::vector<double> beta_vec(beta_obs.begin(), beta_obs.end());
  double beta_mean = 1.0;
  if (adaptive_lambda && beta_vec.size() == D) {
    double sum = 0.0;
    for (int i = 0; i < D; i++) sum += std::abs(beta_vec[i]);
    beta_mean = sum / D;
  }

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (int t = 0; t < nthreads; t++) {
    threads.emplace_back([&, t]() {
      for (int d = t; d < D; d += nthreads) {
        // Row-wise adaptive penalty: lambda_p depends on beta_obs[p], shared across columns
        VectorXd lambda_vec(D);
        for (int p = 0; p < D; ++p) {
          double lambda_p = lambda;
          if (adaptive_lambda && beta_vec.size() == D) {
            lambda_p *= std::abs(beta_vec[p]) / (beta_mean + 1e-8);
          }
          lambda_vec[p] = lambda_p / rho;
        }
        V_next.col(d) = lasso_weighted(U, glm_Y.col(d), V.col(d), lambda_vec, solve_its, fixd ? d + 1 : -1);
      }
    });
  }
  for (auto& thread : threads) thread.join();
  return V_next;
}
