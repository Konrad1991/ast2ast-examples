dot_product <- function(v1, v2) {
  result = 0.0;
  for (i in 1:length(v1)) {
    result = result + v1[i] * v2[i];
  }
  return(result)
}

dot_product <- function(v1, v2) {
  result = 0.0;
  sz <- length(v1)
  i <- 1
  while(i <= sz) {
    result = result + v1[i] * v2[i];
    i <- i + 1
  }
  return(result)
}

v1 <- seq(1, 10^3, 1)
v2 <- runif(n = 10^3)

dot_product_cpp <- function(v1, v2) {
  result::double = 0.0;
  sz::int <- length(v1)
  i::int <- 1
  while(i <= sz) {
    result = result + at(v1, i) * at(v2, i);
    i <- i + 1
  }
  return(result)
}

dp_cpp <- ast2ast::translate(dot_product_cpp,
  handle_inputs = "borrow", output = "R",
  verbose = TRUE
)
cat("\n")

dot_product(v1, v2)
dp_cpp(v1, v2)
microbenchmark::microbenchmark(
  dot_product(v1, v2),
  dp_cpp(v1, v2)
)
stop("bla")





setwd("/home/konrad/Documents/GitHub/RProjects/ast2ast-examples/ParticleSwarmOptimization")
Rcpp::sourceCpp("PSO.cpp")

rosenbrock <- function(parameter) {
  value <- 0
  for (i in 1:(length(parameter) - 1)) {
    value <- value +
      100 * (parameter[i + 1] - parameter[i]^2)^2 +
      (1 - parameter[i])^2
  }
  return(value)
}


lb <- -10000
ub <- 10000
error_threshold <- 0.0000001
npop <- 40

pso(rep(lb, 3), rep(ub, 3), rosenbrock, 10000, npop, error_threshold)


rosenbrock_cpp <- function(parameter) {
  value::double <- 0
  sz::int <- length(parameter)
  for (i in 1:(sz - 1)) {
    value <- value +
      100 * (parameter[i + 1] - parameter[i]^2)^2 +
      (1 - parameter[i])^2
  }
  return(value)
}

rosenbrock_cpp <- ast2ast::translate(rosenbrock_cpp,
  output = "XPtr",
  handle_inputs = "",
  references = TRUE,
  verbose = TRUE
)

pso_xptr(rep(lb, 3), rep(ub, 3),
  rosenbrock_cpp, 1000,
  npop, error_threshold
)
pso_xptr_parallel(rep(lb, 3), rep(ub, 3),
  rosenbrock_cpp, 1000,
  npop, error_threshold
)

r_fct <- function() {
  set.seed(1234)
  pso(rep(lb, 3), rep(ub, 3), rosenbrock, 1000, npop, error_threshold)
}

cpp_fct <- function() {
  set.seed(1234)
  pso_xptr(rep(lb, 3), rep(ub, 3), rosenbrock_cpp, 1000, npop, error_threshold)
}

cpp_fct_parallel <- function() {
  set.seed(1234)
  pso_xptr_parallel(rep(lb, 3), rep(ub, 3), rosenbrock_cpp, 1000, npop, error_threshold)
}


res <- microbenchmark::microbenchmark(
  r_fct(),
  cpp_fct(),
  cpp_fct_parallel()
)
res
getwd()
pdf("RosenbrockBenchmark.pdf")
boxplot(res)
dev.off()



