lambda_fct <- function(u, c) {
  return (5000 * exp(-(10 * u - 3) ^ 2 / 2) + c)
}

lambda_vec_1 <- lambda_fct((1:t_len) / t_len, 1000)

lambda_fct2 <- function(u, c) {
  return (6000 * exp(-(10 * u - 3) ^ 2 / 2) + c)
}

lambda_vec_2 <- lambda_fct2((1:t_len) / t_len, 1000)

plot(1:t_len, lambda_vec_1,  ylim = c(0, max(lambda_vec_1, lambda_vec_2) + 100), xlab="u", ylab = "", type = "l")
lines(1:t_len, lambda_vec_2, type = "l", col = "red")
