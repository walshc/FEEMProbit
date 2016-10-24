\name{FEEMProbit}
\alias{FEEMProbit}
\title{Estimate the individual fixed effects expectation-maximization probit
  model in Chen (2016)}
\description{Function to estimate the individual fixed effects
  expecation-maximization estimator in Chen (2016)}
\usage{
FEEMProbit(formula, data, tol = 1e-9, show.progress = FALSE)
}
\arguments{
  \item{formula}{A symbolic description of the model to be estimated. The
    formula has three parts, the dependent variable, the independent variable
    and the fixed effect. For example, \code{y ~ x1 + x2 | id}, where \code{id}
    is the fixed effect}
  \item{data}{A \code{data.frame}.}
  \item{tol}{Tolerance level to use for the EM algorithm. Lower tolerances
    increase computation time.}
  \item{show.progress}{If \code{TRUE}, print the norm between parameter vectors
    in between each iteration.}
}
\value{
  \item{call}{The function call.}
  \item{coefficients}{The estimated coefficients.}
  \item{fixed.effects}{The estimated individual fixed effects.}
  \item{fitted.values}{The fitted values of the model.}
  \item{residuals}{The model's residuals.}
  \item{model}{The \code{data.frame} used to estimate the model.}
}
\examples{
library(data.table)
set.seed(1)

# Create data:
N <- 200
nT <- 50
df <- data.table(id = rep(1:N, each = nT),
                 t  = rep(1:nT, N))

# Individual fixed effects:
alpha <- (1:N - N/2)/(N/2)
df[, alpha := rep(alpha, each = nT)]

# Explanatory variable:
df[, x := rnorm(N * nT, mean = -2, sd = 0.5)]

# Coefficient to be estimated:
beta <- 0.5

# Latent variable:
df[, ystar := beta * x + alpha + rnorm(N * nT)]

# Observed variable:
df[, y := as.numeric(ystar > 0)]

# Estimate the model:
FEEMProbit(y ~ x | id, data = df, tol = 1e-6)
}
