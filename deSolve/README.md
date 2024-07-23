## Using ast2ast in deSolve

- experimenting with ast2ast in deSolve

```R
#Rcpp::compileAttributes("dsa2a")
install.packages("dsa2a", type = "source", repos = NULL)
.rs.restartR()

library(dsa2a)
library(deSolve)
library(scatterplot3d)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X + Y * Z
    dY <- b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

parameters <- c(a = -8/3, b = -10, c =  28)
state <- c(X = 1, Y = 1, Z = 1)
times <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
#plot(out)
#scatterplot3d(out[,-1], type="l")

out <- ode(state, times, func = "derivs", parms = parameters,
           dllname = "dsa2a", initfunc = "initmod")
#plot(out)
#scatterplot3d(out[,-1], type="l")

out <- ode(state, times, func = "derivs_a2a", parms = parameters,
           dllname = "dsa2a", initfunc = "initmod_a2a")
#plot(out)
#scatterplot3d(out[,-1], type="l")

microbenchmark::microbenchmark(
  out <- ode(y = state, times = times, func = Lorenz, parms = parameters),
  out <- ode(state, times, func = "derivs", parms = parameters, dllname = "dsa2a", initfunc = "initmod"),
  out <- ode(state, times, func = "derivs_a2a", parms = parameters, dllname = "dsa2a", initfunc = "initmod_a2a"))
```
