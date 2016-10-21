FEEMProbit <- function(formula, data, tol = 1e-9, show.progress = FALSE) {

  cl <- match.call()

  # Check if options were used correctly:
  if (!("data.frame" %in% class(data))) {
    stop("'data' must be a 'data.frame'")
  }

  if (!(is.numeric(tol))) {
    stop("'tol' must be numeric. Default is 1e-9.")
  }

  mf <- model.frame(as.Formula(formula), data)
  mf$orig.row.id <- 1:nrow(mf)
  mf <- mf[order(df[[ncol(mf) - 1]]), ]
  id <- as.matrix(mf[, ncol(mf) - 1])
  y <- as.matrix(mf[[1]])
  X <- as.matrix(mf[, c(2:(ncol(mf) - 2))])
  N <- length(unique(id))

  times <- aggregate(y ~ id, FUN = length)$V1

  # Functions to use inside E-M algorithm:
  get_mu_k <- function(X, beta_k, alpha_k) {
    X %*% beta_k + rep(alpha_k, times = times)
  }
  get_y_k <- function(y, mu_k) {
    p <- pnorm(mu_k)
    mu_k + ((y - p) * dnorm(mu_k)) / (p * (1 - p))
  }
  get_beta_k <- function(y_k, alpha_k) {
    solve(crossprod(X)) %*% crossprod(X, y_k - rep(alpha_k, times = times))
  }
  get_alpha_k <- function(y_k, beta_k) {
    df <- data.table(a = c(y_k - X %*% beta_k), id = rep(1:N, times = times))
    df[, mean(a), by = id]$V1
  }

  # Run the E-M algorithm:
  dist <- tol + 1
  beta_k <- matrix(rep(0, ncol(X)))
  alpha_k <- matrix(rep(0, N))
  while (dist > tol) {
    # E-step:
    mu_k <- get_mu_k(X, beta_k, alpha_k)
    y_k <- get_y_k(y, mu_k)
    # M-step:
    new_beta_k <- get_beta_k(y_k, alpha_k)
    new_alpha_k <- get_alpha_k(y_k, new_beta_k)
    dist <- norm(as.matrix(c(new_beta_k, new_alpha_k) - c(beta_k, new_alpha_k)))
    beta_k <- new_beta_k
    alpha_k <- new_alpha_k
    if (show.progress) {
      print(dist)
    }
  }
  beta_k <- c(beta_k)
  names(beta_k) <- names(mf)[2:(ncol(mf)-2)]
  names(alpha_k) <- unique(id)
  mf$predicted.values <- c(X %*% beta_k + rep(alpha_k, times = times))
  mf <- mf[order(mf$orig.row.id), ]
  mf$orig.row.id <- NULL
  result <- list(call = cl, coefficients = beta_k)
  result$fixed.effects <- data.frame(id = sort(unique(id)), fe = alpha_k)
  result$predicted.values <- mf$predicted.values
  result$model <- mf
  class(result) <- "FEEMProbit"
  return(result)
}

print.FEEMProbit <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
