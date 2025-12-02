#' Simula un proceso VECM(0) según la ecuación 3 del artículo Bykhovskaya & Gorin (2022)
#'
#' @param T Número de observaciones temporales
#' @param N Número de series (por defecto 10)
#' @return Una lista con:
#'   - X: matriz de datos simulados (T x N)
#'   - beta_true: vector beta verdadero
#'   - alpha_true: vector alpha verdadero
#'   - Pi: matriz Pi = alpha * beta'
#' @export
simulate_vecm_ec3 <- function(T, N = 10) {
  # Parámetros verdaderos - solo las primeras 2 series están cointegradas
  beta_true <- c(1, -1, rep(0, N-2))
  alpha_true <- c(-0.5, 0.3, rep(0, N-2))
  
  # Matriz Pi = alpha * beta' según el modelo VECM(0)
  Pi <- alpha_true %*% t(beta_true)
  
  # Simular series según la ecuación: ΔX_t = Pi * X_{t-1} + ε_t
  X <- matrix(0, nrow = T, ncol = N)
  X[1, ] <- rnorm(N, 0, 1)
  epsilon <- matrix(rnorm(T * N, 0, 0.2), nrow = T, ncol = N)
  
  for (t in 2:T) {
    X[t, ] <- X[t-1, ] + Pi %*% X[t-1, ] + epsilon[t, ]
  }
  
  return(list(X = X, beta_true = beta_true, alpha_true = alpha_true, Pi = Pi))
}
