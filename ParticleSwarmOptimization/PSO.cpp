// TODO: add parallel pso in order to show that translated functions can be called in parallel

// [[Rcpp::depends(ast2ast)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::plugins(cpp2a)]]
#include "etr.hpp"
#include <RcppThread.h>

using namespace Rcpp;
using namespace arma;

// Define the PSO function
// [[Rcpp::export]]
List pso(const vec &lb, const vec &ub, Function loss, int ngen, int npop,
         double error_threshold) {
  int npar = lb.n_elem;
  mat swarm = randu<mat>(npop, npar);
  mat v = zeros<mat>(npop, npar);
  vec swarm_bests = zeros<vec>(npop);
  vec swarm_errors = zeros<vec>(npop);

  const double initial_cog = 2.5;
  const double final_cog = 0.5;
  const double initial_soc = 0.5;
  const double final_soc = 2.5;
  const double w_max = 0.9;
  const double w_min = 0.4;

  for (int i = 0; i < npop; ++i) {
    swarm.row(i) = (randu<rowvec>(npar) % (ub.t() - lb.t())) + lb.t();
    swarm_errors[i] = as<double>(loss(swarm.row(i).t()));
    swarm_bests[i] = swarm_errors[i];
  }

  uword global_best = index_min(swarm_bests);
  rowvec global_best_vec = swarm.row(global_best);
  double global_best_error = swarm_bests[global_best];
  mat swarm_best_params = swarm;

  auto correctBounds = [](const rowvec &target, const vec &lb, const vec &ub) {
    rowvec res(target.size());
    for (size_t i = 0; i < lb.size(); i++) {
      if (target(i) < lb(i)) {
        res(i) = lb(i);
      } else if (target(i) > ub(i)) {
        res(i) = ub(i);
      } else {
        res(i) = target(i);
      }
    }
    return res;
  };

  int iter = 0;
  while (iter < ngen) {
    double w = w_max - iter * (w_max - w_min) / ngen;
    double cog = initial_cog - (initial_cog - final_cog) * (iter + 1) / ngen;
    double soc = initial_soc - (initial_soc - final_soc) * (iter + 1) / ngen;

    for (int i = 0; i < npop; ++i) {
      rowvec local_best_vec = swarm.row(global_best);

      v.row(i) =
          w * v.row(i) +
          cog * randu<double>() * (swarm_best_params.row(i) - swarm.row(i)) +
          soc * randu<double>() * (local_best_vec - swarm.row(i));
      swarm.row(i) += v.row(i);

      swarm.row(i) = correctBounds(swarm.row(i), lb, ub);

      double error = as<double>(loss(swarm.row(i).t()));

      if (!std::isinf(error) && !std::isnan(error) && error < swarm_bests[i]) {
        swarm_bests[i] = error;
        swarm_best_params.row(i) = swarm.row(i);
      }
      if (!std::isinf(error) && !std::isnan(error) &&
          error < global_best_error) {
        global_best = i;
        global_best_vec = swarm.row(i);
        global_best_error = error;
      }
    }

    iter++;

    if (global_best_error < error_threshold) {
      break;
    }
  }

  double insilico = as<double>(loss(global_best_vec.t()));
  return List::create(Named("insilico") = insilico,
                      Named("global_best_vec") = global_best_vec);
}

typedef double (*FP)(etr::Vec<double> &a);

// [[Rcpp::export]]
List pso_xptr(const vec &lb, const vec &ub, Rcpp::XPtr<FP> loss_xptr, int ngen,
              int npop, double error_threshold) {
  FP loss = *loss_xptr;
  int npar = lb.n_elem;
  mat swarm = randu<mat>(npop, npar);
  mat v = zeros<mat>(npop, npar);
  vec swarm_bests = zeros<vec>(npop);
  vec swarm_errors = zeros<vec>(npop);

  const double initial_cog = 2.5;
  const double final_cog = 0.5;
  const double initial_soc = 0.5;
  const double final_soc = 2.5;
  const double w_max = 0.9;
  const double w_min = 0.4;

  etr::Vec<double> param(etr::SI{npar});

  for (int i = 0; i < npop; ++i) {
    swarm.row(i) = (randu<rowvec>(npar) % (ub.t() - lb.t())) + lb.t();
    for (size_t j = 0; j < npar; j++) {
      param[j] = swarm.row(i)(j);
    }
    swarm_errors[i] = loss(param);
    swarm_bests[i] = swarm_errors[i];
  }

  uword global_best = index_min(swarm_bests);
  rowvec global_best_vec = swarm.row(global_best);
  double global_best_error = swarm_bests[global_best];
  mat swarm_best_params = swarm;

  auto correctBounds = [](const rowvec &target, const vec &lb, const vec &ub) {
    rowvec res(target.size());
    for (size_t i = 0; i < lb.size(); i++) {
      if (target(i) < lb(i)) {
        res(i) = lb(i);
      } else if (target(i) > ub(i)) {
        res(i) = ub(i);
      } else {
        res(i) = target(i);
      }
    }
    return res;
  };

  int iter = 0;
  while (iter < ngen) {
    double w = w_max - iter * (w_max - w_min) / ngen;
    double cog = initial_cog - (initial_cog - final_cog) * (iter + 1) / ngen;
    double soc = initial_soc - (initial_soc - final_soc) * (iter + 1) / ngen;

    for (int i = 0; i < npop; ++i) {
      rowvec local_best_vec = swarm.row(global_best);

      v.row(i) =
          w * v.row(i) +
          cog * randu<double>() * (swarm_best_params.row(i) - swarm.row(i)) +
          soc * randu<double>() * (local_best_vec - swarm.row(i));
      swarm.row(i) += v.row(i);

      swarm.row(i) = correctBounds(swarm.row(i), lb, ub);
      for (size_t j = 0; j < npar; j++) {
        param[j] = swarm.row(i)(j);
      }

      double error = loss(param);

      if (!std::isinf(error) && !std::isnan(error) && error < swarm_bests[i]) {
        swarm_bests[i] = error;
        swarm_best_params.row(i) = swarm.row(i);
      }
      if (!std::isinf(error) && !std::isnan(error) &&
          error < global_best_error) {
        global_best = i;
        global_best_vec = swarm.row(i);
        global_best_error = error;
      }
    }

    iter++;

    if (global_best_error < error_threshold) {
      break;
    }
  }

  for (size_t j = 0; j < npar; j++) {
    param[j] = global_best_vec(j);
  }
  double insilico = loss(param);
  return List::create(Named("insilico") = insilico,
                      Named("global_best_vec") = global_best_vec);
}





// [[Rcpp::export]]
List pso_xptr_parallel(const vec &lb, const vec &ub, Rcpp::XPtr<FP> loss_xptr, int ngen,
              int npop, double error_threshold) {
  FP loss = *loss_xptr;
  int npar = lb.n_elem;
  mat swarm = randu<mat>(npop, npar);
  mat v = zeros<mat>(npop, npar);
  vec swarm_bests = zeros<vec>(npop);
  vec swarm_errors = zeros<vec>(npop);

  const double initial_cog = 2.5;
  const double final_cog = 0.5;
  const double initial_soc = 0.5;
  const double final_soc = 2.5;
  const double w_max = 0.9;
  const double w_min = 0.4;

  std::vector<etr::Vec<double>> params(npop);
  for (size_t i = 0; i < params.size(); i++) {
    params[i] = etr::Vec<double>(etr::SI{npar});
  }

  for (int i = 0; i < npop; ++i) {
    swarm.row(i) = (randu<rowvec>(npar) % (ub.t() - lb.t())) + lb.t();
    for (size_t j = 0; j < npar; j++) {
      params[i][j] = swarm.row(i)(j);
    }
    swarm_errors[i] = loss(params[i]);
    swarm_bests[i] = swarm_errors[i];
  }

  uword global_best = index_min(swarm_bests);
  rowvec global_best_vec = swarm.row(global_best);
  double global_best_error = swarm_bests[global_best];
  mat swarm_best_params = swarm;

  int number_threads = std::thread::hardware_concurrency() - 1;
  RcppThread::ThreadPool pool(number_threads);
  std::vector<std::future<double>> futures(npop);

  auto correctBounds = [](const rowvec &target, const vec &lb, const vec &ub) {
    rowvec res(target.size());
    for (size_t i = 0; i < lb.size(); i++) {
      if (target(i) < lb(i)) {
        res(i) = lb(i);
      } else if (target(i) > ub(i)) {
        res(i) = ub(i);
      } else {
        res(i) = target(i);
      }
    }
    return res;
  };

  int iter = 0;
  while (iter < ngen) {
    double w = w_max - iter * (w_max - w_min) / ngen;
    double cog = initial_cog - (initial_cog - final_cog) * (iter + 1) / ngen;
    double soc = initial_soc - (initial_soc - final_soc) * (iter + 1) / ngen;

    for (int i = 0; i < npop; ++i) {
      rowvec local_best_vec = swarm.row(global_best);

      v.row(i) =
          w * v.row(i) +
          cog * randu<double>() * (swarm_best_params.row(i) - swarm.row(i)) +
          soc * randu<double>() * (local_best_vec - swarm.row(i));
      swarm.row(i) += v.row(i);

      swarm.row(i) = correctBounds(swarm.row(i), lb, ub);
    }

    for (int i = 0; i < npop; ++i) {
      for (size_t j = 0; j < npar; j++) {
        params[i][j] = swarm.row(i)(j);
      }
    }

    for (int i = 0; i < npop; ++i) {
      futures[i] = pool.pushReturn(loss, std::ref(params[i]));
    }

    for (int i = 0; i < npop; ++i) {
      swarm_errors[i] = futures[i].get();
    }

    for (int i = 0; i < npop; i++) {
      if (!std::isinf(swarm_errors[i]) &&
        !std::isnan(swarm_errors[i]) &&
        swarm_errors[i]< swarm_bests[i]) {
        swarm_bests[i] = swarm_errors[i];
        swarm_best_params.row(i) = swarm.row(i);
      }
      if (!std::isinf(swarm_errors[i]) && !std::isnan(swarm_errors[i]) &&
          swarm_errors[i] < global_best_error) {
        global_best = i;
        global_best_vec = swarm.row(i);
        global_best_error =swarm_errors[i];
      }
    }

    iter++;

    if (global_best_error < error_threshold) {
      break;
    }
  }

  for (size_t j = 0; j < npar; j++) {
    params[0][j] = global_best_vec(j);
  }
  double insilico = loss(params[0]);
  return List::create(Named("insilico") = insilico,
                      Named("global_best_vec") = global_best_vec);
}

