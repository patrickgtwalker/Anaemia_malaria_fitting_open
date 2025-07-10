# Function to compute cubic spline interpolation
cubic_spline <- function(x, y, x_pred) {
  # Number of intervals
  n <- length(x) - 1
  
  # Initialize arrays
  c <- l <- mu <- z <- rep(0, n + 1)
  h <- b <- d <- alpha <- rep(NA, n)
  
  # Compute the differences between consecutive x values
  for (i in 1:n) {
    h[i] <- x[i + 1] - x[i]
  }
  
  # Compute alpha values for the system of equations
  for (i in 2:n) {
    alpha[i] <- 3/h[i] * (y[i + 1] - y[i]) - 3/h[i - 1] * (y[i] - y[i - 1])
  }
  
  # Forward sweep to solve for c values
  l[1] <- 1
  for (i in 2:n) {
    l[i] <- 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
    mu[i] <- h[i]/l[i]
    z[i] <- (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
  }
  l[n + 1] <- 1
  
  # Backward substitution to calculate the coefficients
  for (i in n:1) {
    c[i] <- z[i] - mu[i] * c[i + 1]
    b[i] <- (y[i + 1] - y[i])/h[i] - h[i] * (c[i + 1] + 2 * c[i])/3
    d[i] <- (c[i + 1] - c[i])/(3 * h[i])
  }
  
  # Calculate the spline for each prediction point
  s <- rep(NA, length(x_pred))
  j <- 1
  for (i in seq_along(x_pred)) {
    while (x_pred[i] > x[j + 1]) {
      j <- j + 1
    }
    s[i] <- y[j] + b[j] * (x_pred[i] - x[j]) + c[j] * (x_pred[i] - x[j])^2 + d[j] * (x_pred[i] - x[j])^3
  }
  
  return(s)
}
