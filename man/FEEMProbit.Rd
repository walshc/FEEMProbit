\name{FEEMProbit}
\alias{FEEMProbit}
\title{Estimate the individual fixed effects expectation-maximization probit
  model in Chen (2016)}
\description{Function to estimate the individual fixed effects expecation-maximization estimator in Chen (2016)}
\usage{
FEEMProbit(formula, data, tol = 1e-9, show.progress = FALSE)
}
\arguments{
  \item{formula}{A symbolic description of the model to be estimated.}
  \item{data}{A \code{pdata.frame}.}
  \item{tol}{Tolerance level to use for the EM algorithm. Lower tolerances
    increase computation time.}
  \item{show.progress}{If \code{TRUE}, print the norm between parameter vectors
    in between each iteration.}
}
\details{
  A \code{pdata.frame} from the package \code{plm} must be used. The fixed
  effects that will be estimated are the \code{id} variables in that
  \code{pdata.frame}. They do not need to be specified in the formula.
}
\value{
  \item{call}{The function call.}
  \item{coefficients}{The estimated coefficients.}
  \item{fixed.effects}{The estimated individual fixed effects.}
}
\examples{
library(data.table)
library(plm)
set.seed(1)
N <- 100
nT <- 50

df <- data.table(id = rep(1:N, each = nT),
                 t  = rep(1:nT, N))

# Individual fixed effects:
alpha <- rnorm(N, sd = 0.2)
df[, alpha := rep(alpha, each = nT)]

# Explanatory variable:
df[, x := rnorm(N * nT, mean = -2, sd = 0.5)]

# Coefficient to be estimated:
beta <- 0.5

# Latent variable:
df[, ystar := x * beta + alpha + rnorm(N * nT)]

# Observed variable:
df[, y := as.numeric(ystar > 0)]

# Convert to pdata.frame:
df <- pdata.frame(df, index = c("id", "t"))

# Estimate the model:
FEEMProbit(y ~ x, data = df)
}
