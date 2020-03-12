#' @title Teste da razão de verossimilhança generalizada
#' @author Pedro Rafael D. Marinho
#' @description Calcula a estatística da razão de verossimilhança generalizada.
#' @details A função recebe como argumento uma função densidade de probabilidade ou uma função de probabilidade
#' que é passada como argumento à \code{f}. O objeto passado como argumento de \code{f} deverá ser implementado conforme os
#' exemplos. Note que se a função for implementada segundo os exemplos, não haverá a necessidade de implementar a função sob a
#' hipótese nula, visto que são as mesmas funções. Para especificar os parâmetros que serão fixados, i.e, para especificar
#' a distribuição sob a hipótese nula, utiliza-se o argumento \code{par0} que receberá uma lista formada por dois vetores.
#' O primeiro vetor da lista de verá ser um vetor de strings com os nomes das variáveis que deseja-se fixar e o segundo vetor
#' deverá conter os valores que serão atribuidos à cada uma das variáveis passada ao primeiro vetor.
#' @param f Função densidade de probabilidade considerada no teste da razão de verossimilhança. Essa função deverá ser
#' implementada conforme os exemplos.
#' @param data Um conjunto de dados que será considerado para a realização do teste.
#' @param par0 Uma lista contendo como primeiro elemento um vetor com o nomes das variáveis que serão fixadas como hipótese nula
#' (variáveis aos quais desejamos testar) e um vetor com os valores fixados para cada uma das respectivas variáveis.
#' @param par1 Lista com dois elementos, sendo o primeiro um vetor com os nomes das variáveis que receberão valores fixos sob a
#' hipótese alternativa e o segundo elemendo é um outro vetor com os valores impostos às variáveis.
#' @param kicks Chutes iniciais que desejamos considerar no método de otimização.
#' @param ... Lista de argumentos adicionais que serão passados à função \code{optim()} otimizada para otimização. Por exemplo,
#' será possível escolher o método de otimização a ser utilizado.
#' @importFrom stats optim
#' @importFrom stats qchisq
#' @importFrom magrittr "%>%"
#' @return Valor da estatística do teste da razão de verossimilhança generalizada.
#' @examples
#' pdf_w <- function(par, x, var = NULL){
#'    alpha <- par[1]
#'    beta <- par[2]
#'    dweibull(x, shape = alpha, scale = beta)
#' }
#'
#' rw <- function(n = 1L, alpha, beta){
#'    rweibull(n = n, shape = alpha, scale = beta)
#' }
#'
#' data <- rw(n = 100L, alpha = 1, beta = 1)
#'
#' lrt(f = pdf_w, data = data, kicks = c(1, 1), par0 = list("beta", 1))
#' @export
lrt <- function(f, data, kicks, par0 = NULL, ...) {
  if (is.null(par0))
    stop("Informar uma lista informando o parâmetro e o valor sob a hipótese nula.")

  body(f) %<>% as.list %>%
    append(quote(if (is.list(var))
      eval(parse(
        text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")
      ))), length(body(f)) - 1L) %>%
    as.call %>%
    as.expression

  # Log-Likelihood under the null hypothesis. -----------------------------------
  log_lik_h0 <- function(par, x) {
    -sum(log(f(par, x, var = par0)))
  }

  # Unrestricted log-likelihood. -----------------------------------
  log_lik <- function(par, x) {
    -sum(log(f(par, x)))
  }

  myoptim <-
    function(...)
      tryCatch(
        expr = optim(...),
        error = function(e)
          NA
      )

  par_h0 <- myoptim(par = kicks, fn = log_lik_h0, x = data, ...)

  if (!is.list(par_h0) || par_h0$convergence != 0L)
    return(NA)

  par_h <- myoptim(par = kicks, fn = log_lik, x = data, ...)

  if (!is.list(par_h) || par_h$convergence != 0L)
    return(NA)

  lambda <-
    2 * (log_lik_h0(par = par_h0$par, x = data) - log_lik(par_h$par, x = data))

  lambda[lambda < 0] <- 0

  # Estatística de razão de verossimilhança:
  lambda
}
