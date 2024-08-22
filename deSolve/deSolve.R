library(deSolve)
library(scatterplot3d)
is_installed <- function(pkg) {
  is.element(pkg, installed.packages()[, 1])
}

# R solution
# ==========================================================================
Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X + Y * Z
    dY <- b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

parameters <- c(a = -8 / 3, b = -10, c = 28)
state <- c(X = 1, Y = 1, Z = 1)
times <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
plot(out)
scatterplot3d(out[, -1], type = "l")

# ast2ast
# ==========================================================================
lorenz_fct <- function(t, state, state_deriv) {
  state_deriv[1] <- a * state[1] + state[2] * state[3]
  state_deriv[2] <- b * (state[2] - state[3])
  state_deriv[3] <- -state[1] * state[2] + c * state[2] - state[3]
}
odefcpp <- ast2ast::translate(
  f = lorenz_fct,
  types_of_args = rep("double", 3),
  data_structures = c("scalar", "borrow", "borrow"),
  handle_inputs = rep("", 3),
  references = c(FALSE, TRUE, TRUE),
  output = "XPtr",
  getsource = TRUE
)
odefcpp <- strsplit(odefcpp, "\n")[[1]]
odetmp <- character(length(odefcpp))
counter <- 1
fct_start <- which(grepl("void lorenz_fct", odefcpp))
repeat {
  odetmp[counter] <- odefcpp[[fct_start]]
  if (odefcpp[[fct_start]] == "}") break
  fct_start <- fct_start + 1
  counter <- counter + 1
}
odetmp <- paste(odetmp[1:counter], collapse = "\n")

code <- '
// [[Rcpp::depends(ast2ast, RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

#include "etr.hpp"
using namespace Rcpp;
using namespace etr;

static double parms[3];
#define a parms[0]
#define b parms[1]
#define c parms[2]

extern "C" {
void initmod_ast2ast(void (*odeparms)(int *, double *));
}

extern "C" {
void derivs_ast2ast(int *neq, double *t, double *y, double *ydot, double *yout,
                    int *ip);
}

void initmod_ast2ast(void (*odeparms)(int *, double *)) {
  int N = 3;
  odeparms(&N, parms);
}

%s

/* Derivatives */
void derivs_ast2ast(int *neq, double *t, double *y, double *ydot, double *yout,
                    int *ip) {
  BorrowPtr y_(y, *neq);
  BorrowPtr ydot_(ydot, *neq);
  double t_ = *t;
  lorenz_fct(t_, y_, ydot_);
}
'
odetmp <- sprintf(code, odetmp)

Rcpp::sourceCpp(code = odetmp, cacheDir = ".", verbose = TRUE)
dir <- list.dirs()[3] # not the best way
shared_objects <- list.files(dir, pattern = "\\.so$", full.names = TRUE) # dll on windows
dyn.load(shared_objects)

out <- ode(state, times,
  func = "derivs_ast2ast", parms = parameters,
  dllname = "sourceCpp_2", initfunc = "initmod_ast2ast"
)
plot(out)
scatterplot3d(out[, -1], type = "l")

# plain c
# ==========================================================================
system("R CMD SHLIB 'pure.c'")
dyn.load("pure.so")
out <- ode(state, times,
  func = "derivs", parms = parameters,
  dllname = "pure", initfunc = "initmod"
)
plot(out)
scatterplot3d(out[, -1], type = "l")

microbenchmark::microbenchmark(
  out <- ode(
    y = state,
    times = times,
    func = Lorenz,
    parms = parameters
  ),
  out <- ode(state, times,
    func = "derivs_ast2ast",
    parms = parameters,
    dllname = "sourceCpp_2",
    initfunc = "initmod_ast2ast"
  ),
  out <- ode(state, times,
    func = "derivs",
    parms = parameters,
    dllname = "pure",
    initfunc = "initmod"
  )
)
