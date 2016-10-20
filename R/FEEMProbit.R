FEEMProbit <- function(formula, data, tol = 1e-9, show.progress = FALSE) {

  cl <- match.call()

  # Check if options were used correctly:
  if (!("pdata.frame" %in% class(data))) {
    stop("'data' must be a 'pdata.frame'. Use package 'plm'.")
  }

  if (!(is.numeric(tol))) {
    stop("'tol' must be numeric. Default is 1e-9.")
  }

  # Check if data is balanced and if not balance it:
  if (!plm::pdim(data)$balanced) {
    un.id <- sort(unique(index(data, "id")))
    un.time <- sort(unique(index(data, "time")))
    rownames(data) <- paste(index(data, "id"), index(data, "time"), sep = ".")
    allRows <- as.character(t(outer(un.id, un.time, paste, sep = ".")))
    data <- data[allRows, ]
    rownames(data) <- allRows
    index <- data.frame(id = rep(un.id, each = length(un.time)),
                        time = rep(un.time, length(un.id)),
                        row.names = rownames(data))
    class(index) <- c("pindex", "data.frame")
    attr(data, "index") <- index
    msg <- paste("This function doesn't work with unbalanced panel data.",
                 "Forcing to be balanced")
    warning(msg)
  }

  N <- length(unique(index(data, "id")))
  nT <- length(unique(index(data, "time")))
  mf <- model.frame(formula, data)
  y <- as.matrix(mf[[1]])
  X <- as.matrix(mf[, c(2:ncol(mf))])

  # Functions to use inside E-M algorithm:
  get_mu_k <- function(X, beta_k, alpha_k) {
    X %*% beta_k + rep(alpha_k, each = nT)
  }
  get_y_k <- function(y, mu_k) {
    p <- pnorm(mu_k)
    mu_k + ((y - p) * dnorm(mu_k)) / (p * (1 - p))
  }
  get_beta_k <- function(y_k, alpha_k) {
    solve(crossprod(X)) %*% crossprod(X, y_k - rep(alpha_k, each = nT))
  }
  get_alpha_k <- function(y_k, beta_k) {
    df <- data.table(a = c(y_k - X %*% beta_k), id = rep(1:N, each = nT))
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
  names(beta_k) <- names(mf)[2:ncol(mf)]
  names(alpha_k) <- sort(unique(index(data, "id")))
  predicted.values <- c(X %*% beta_k + rep(alpha_k, each = nT))
  result <- list(call = cl, coefficients = beta_k, fixed.effects = alpha_k,
                 predicted.values = predicted.values)
  class(result) <- "FEEMProbit"
  return(result)
}

print.FEEMProbit <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
