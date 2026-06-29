# MLG com Distribuição Bell ----

## Pacotes
# Distribuição Bell
require(bellreg)
# MLG Binomial Negativa
require(MASS)
# Adicionais MLG
require(car)
# Visualização
require(knitr)
require(ggplot2)
require(RColorBrewer)


## Funções e Variáveis Globais

# Funções
## resíduos quantílicos escalonados do MLG Bell
scaled_residuals_bell = function(out, random = TRUE, B = 250){
  # out = objeto glm com distribuição Bell
  # random = aleatorizar ou não os resíduos
  # B = número de replicações bootstrap
  
  ## variável resposta e tamanho da amostra
  vY = model.frame(out)[[1]]
  n = length(vY)
  
  ## predito do modelo
  vMu_hat = predict(out, type = 'response')
  
  ## resíduos PIT
  eps = rep(NA, n)
  # loop principal
  for(i in 1:n){
    aux = rbell(B, theta = LambertW::W(vMu_hat[i]))
    E = ecdf(aux)
    a = E(vY[i] - 1)
    b = E(vY[i])
    # sem aleatorização
    if(random == FALSE){
      eps[i] = b
    }
    # com aleatorização
    if(random == TRUE){
      eps[i] = runif(1, min = a, max = b)
    }
  }
  
  ## saída da função
  L = list('scaled_residuals' = eps)
  return(L)
}

## resíduos quantílicos escalonados do MLG Poisson
scaled_residuals_pois = function(out, random = TRUE, B = 250){
  # out = objeto glm com distribuição Poisson
  # random = aleatorizar ou não os resíduos
  # B = número de replicações bootstrap
  
  ## variável resposta e tamanho da amostra
  vY = model.frame(out)[[1]]
  n = length(vY)
  
  ## predito do modelo
  vMu_hat = predict(out, type = 'response')
  
  ## resíduos PIT
  eps = rep(NA, n)
  # loop principal
  for(i in 1:n){
    aux = rpois(B, lambda = vMu_hat[i])
    E = ecdf(aux)
    a = E(vY[i] - 1)
    b = E(vY[i])
    # sem aleatorização
    if(random == FALSE){
      eps[i] = b
    }
    # com aleatorização
    if(random == TRUE){
      eps[i] = runif(1, min = a, max = b)
    }
  }
  
  ## saída da função
  L = list('scaled_residuals' = eps)
  return(L)
}

## resíduos quantílicos escalonados do MLG NB
scaled_residuals_nb = function(out, random = TRUE, B = 250){
  # out = objeto glm com distribuição Binomial Negativa
  # random = aleatorizar ou não os resíduos
  # B = número de replicações bootstrap
  
  ## variável resposta e tamanho da amostra
  vY = model.frame(out)[[1]]
  n = length(vY)
  
  ## predito do modelo
  vMu_hat = predict(out, type = 'response')
  
  ## parâmetro de sobredispersão da NB
  theta = out$theta
  
  ## resíduos PIT
  eps = rep(NA, n)
  # loop principal
  for(i in 1:n){
    aux = MASS::rnegbin(B, mu = vMu_hat[i], theta = theta)
    E = ecdf(aux)
    a = E(vY[i] - 1)
    b = E(vY[i])
    # sem aleatorização
    if(random == FALSE){
      eps[i] = b
    }
    # com aleatorização
    if(random == TRUE){
      eps[i] = runif(1, min = a, max = b)
    }
  }
  
  ## saída da função
  L = list('scaled_residuals' = eps)
  return(L)
}

# Variáveis
H <- 10




# Densidade da Distribuição Bell ----
bellpmf <- function(alphas, max_x) {
  #alphas é um vetor dos valores de alpha para comparação
  #max_x é um número inteiro positivo de limite para o plot
  y <- 0:max_x
  n <- length(alphas)
  pontos <- 16:(16+n-1)
  cores <- brewer.pal(n, "Dark2")
  ylim <- max(dbell(0:max_x, min(alphas)))
  
  plot(y, dbell(y, theta = alphas[1]),
       type = "b",
       pch = pontos[1],
       col = cores[1],
       lty = 2,
       ylim = c(0, ylim),
       xlab = "y",
       ylab = expression(Pr(Y == y)),
       las = 1)
  
  for (i in 2:n) {
    lines(y, dbell(y, theta = alphas[i]),
          type = "b",
          pch = pontos[i],
          col = cores[i],
          lty = 2)
  }
  
  legend("topright", 
         legend = as.expression(lapply(alphas, function(t) bquote(alpha == .(t)))),
         col = cores,
         pch = pontos,
         lty = 2,
         bty = "n",
         cex = 0.9)
}

par(mfrow = c(1, 1))

alphas <- c(0.4, 0.8, 1.2, 1.6) #Pelo menos 3 valores
max_x <- 15

bellpmf(alphas, max_x)



rm(list = setdiff(ls(), c("scaled_residuals_bell", "scaled_residuals_pois",
                          "scaled_residuals_nb", "H")))
# Dados Simulados ----

## Variável Explicativa Uniforme ----
# Função de dados com variável explicativa uniforme
dadosunif <- function(n) {
  mX = matrix(NA, nrow = n, ncol = 2)
  vBeta = matrix(c(1, 2))
  vY = rep(NA, n)
  vMu = rep(NA, n)
  vX1 = rep(1, n)
  vX2_raw = runif(n, min = 0, max = 1)
  vX2_scaled = scale(vX2_raw, 
                     scale = FALSE, center = TRUE) # padronizando vX2
  for(i in 1:n){
    mX[i, 1] = vX1[i]
    mX[i, 2] = vX2_scaled[i]
    vMu[i] = exp(t(mX[i, ]) %*% vBeta)
    vY[i] = rbell(n = 1, theta = LambertW::W(vMu[i]))
  }
  return(list(vY = vY, mX = mX ))
}

# Vetor de cores usados
cores <- brewer.pal(4, "Blues")

# Bases de dados
set.seed(20)
d1 <- dadosunif(50)
d2 <- dadosunif(100)
d3 <- dadosunif(500)
d4 <- dadosunif(1000)


### Visualização ----
par(mfrow = c(2, 2))
datasets <- list(d1, d2, d3, d4)
n_obs <- c(50, 100, 500, 1000)

for (i in 1:4) {
  hist(datasets[[i]]$vY, 
       prob = TRUE,
       main = paste("Histograma de vY\n(n =", n_obs[i], ")"), 
       col = cores[i],
       xlab = "vY",
       ylab = "Densidade")
  
  lines(density(datasets[[i]]$vY))
}


### Modelos MLG ----
######################################### n = 50
# modelo bell
out.belld1 = glm(vY ~ mX[, 2], data = d1,
                 family = bell(link = log))
# modelo pois
out.poisd1 = glm(vY ~ mX[, 2], data = d1,
                 family = poisson(link = log))
######################################### n = 100
# modelo bell
out.belld2 = glm(vY ~ mX[, 2], data = d2,
                 family = bell(link = log))
# modelo pois
out.poisd2 = glm(vY ~ mX[, 2], data = d2,
                 family = poisson(link = log))
######################################### n = 500
# modelo bell
out.belld3 = glm(vY ~ mX[, 2], data = d3,
                 family = bell(link = log))
# modelo pois
out.poisd3 = glm(vY ~ mX[, 2], data = d3,
                 family = poisson(link = log))
######################################### n = 1000
# modelo bell
out.belld4 = glm(vY ~ mX[, 2], data = d4,
                 family = bell(link = log))
# modelo pois
out.poisd4 = glm(vY ~ mX[, 2], data = d4,
                 family = poisson(link = log))


### AIC ----
tabela_aic <- data.frame(
  Amostra = c("n = 50", "n = 100", "n = 500", "n = 1000"),
  Bell    = c(round(AIC(out.belld1), 2),
              round(AIC(out.belld2), 2),
              round(AIC(out.belld3), 2),
              round(AIC(out.belld4), 2)),
  
  Poisson = c(round(AIC(out.poisd1), 2),
              round(AIC(out.poisd2), 2),
              round(AIC(out.poisd3), 2),
              round(AIC(out.poisd4), 2))
)

kable(tabela_aic,
      col.names = c("Amostra (n)", "Modelo Bell", "Modelo Poisson"),
      align = c("l", "c", "c"),
      row.names = FALSE,
      caption = "Comparação do AIC dos Modelos Bell e Poisson de Amostra")



### Estimativa de Betas ----

# n = 50
summary_belld1 <- summary(out.belld1)
coef_belld1 <- summary_belld1$coefficients[, "Estimate"]

summary_poisd1 <- summary(out.poisd1)
coef_poisd1 <- summary_poisd1$coefficients[, "Estimate"]

# n = 100
summary_belld2 <- summary(out.belld2)
coef_belld2 <- summary_belld2$coefficients[, "Estimate"]

summary_poisd2 <- summary(out.poisd2)
coef_poisd2 <- summary_poisd2$coefficients[, "Estimate"]

# n = 500
summary_belld3 <- summary(out.belld3)
coef_belld3 <- summary_belld3$coefficients[, "Estimate"]

summary_poisd3 <- summary(out.poisd3)
coef_poisd3 <- summary_poisd3$coefficients[, "Estimate"]

# n = 1000
summary_belld4 <- summary(out.belld4)
coef_belld4 <- summary_belld4$coefficients[, "Estimate"]

summary_poisd4 <- summary(out.poisd4)
coef_poisd4 <- summary_poisd4$coefficients[, "Estimate"]

tabela_betas <- data.frame(
  Amostra = c("n = 50", "n = 100", "n = 500", "n = 1000"),
  B0_Bell = c(round(coef_belld1[1], 2),
              round(coef_belld2[1], 2),
              round(coef_belld3[1], 2),
              round(coef_belld4[1], 2)),
  
  B0_Pois = c(round(coef_poisd1[1], 2),
              round(coef_poisd2[1], 2),
              round(coef_poisd3[1], 2),
              round(coef_poisd4[1], 2)),
  
  B1_Bell = c(round(coef_belld1[2], 2),
              round(coef_belld2[2], 2),
              round(coef_belld3[2], 2),
              round(coef_belld4[2], 2)),
  
  B1_Pois = c(round(coef_poisd1[2], 2),
              round(coef_poisd2[2], 2),
              round(coef_poisd3[2], 2),
              round(coef_poisd4[2], 2))
)

# Gerando a tabela de Betas
kable(tabela_betas,
      col.names = c("Amostra (n)", 
                    "Beta 0 (Bell)", "Beta 0 (Poisson)",
                    "Beta 1 (Bell)", "Beta 1 (Poisson)"),
      align = c("l", "c", "c", "c", "c"),
      row.names = FALSE,
      caption = "Estimativas dos Parâmetros Beta 0 e Beta 1 por Amostra")



### Análise de Resíduos ----

#### Resíduos PIT ----
set.seed(20)
## Base de Dados d1
# Bell d1
L_bell_d1 = scaled_residuals_bell(out.belld1, random = TRUE, B = 250)
eps_bell_d1 = L_bell_d1$scaled_residuals
# Poisson d1
L_pois_d1 = scaled_residuals_pois(out.poisd1, random = TRUE, B = 250)
eps_pois_d1 = L_pois_d1$scaled_residuals

## Base de Dados d2
# Bell d2
L_bell_d2 = scaled_residuals_bell(out.belld2, random = TRUE, B = 250)
eps_bell_d2 = L_bell_d2$scaled_residuals
# Poisson d2
L_pois_d2 = scaled_residuals_pois(out.poisd2, random = TRUE, B = 250)
eps_pois_d2 = L_pois_d2$scaled_residuals

## Base de Dados d3
# Bell d3
L_bell_d3 = scaled_residuals_bell(out.belld3, random = TRUE, B = 250)
eps_bell_d3 = L_bell_d3$scaled_residuals
# Poisson d3
L_pois_d3 = scaled_residuals_pois(out.poisd3, random = TRUE, B = 250)
eps_pois_d3 = L_pois_d3$scaled_residuals

## Base de Dados d4
# Bell d4
L_bell_d4 = scaled_residuals_bell(out.belld4, random = TRUE, B = 250)
eps_bell_d4 = L_bell_d4$scaled_residuals
# Poisson d4
L_pois_d4 = scaled_residuals_pois(out.poisd4, random = TRUE, B = 250)
eps_pois_d4 = L_pois_d4$scaled_residuals


##### Histogramas ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
hist(eps_bell_d1, breaks = H, probability = TRUE, main = "PIT Bell (n = 50)", col = cores[1])
abline(h = 1, lty = 2, col = 'red')
# Poisson d1
hist(eps_pois_d1, breaks = H, probability = TRUE, main = "PIT Poisson (n = 50)", col = cores[1])
abline(h = 1, lty = 2, col = 'red')

## Base de Dados d2
# Bell d2
hist(eps_bell_d2, breaks = H, probability = TRUE, main = "PIT Bell (n = 100)", col = cores[2])
abline(h = 1, lty = 2, col = 'red')
# Poisson d2
hist(eps_pois_d2, breaks = H, probability = TRUE, main = "PIT Poisson (n = 100)", col = cores[2])
abline(h = 1, lty = 2, col = 'red')

## Base de Dados d3
# Bell d3
hist(eps_bell_d3, breaks = H, probability = TRUE, main = "PIT Bell (n = 500)", col = cores[3])
abline(h = 1, lty = 2, col = 'red')
# Poisson d3
hist(eps_pois_d3, breaks = H, probability = TRUE, main = "PIT Poisson (n = 500)", col = cores[3])
abline(h = 1, lty = 2, col = 'red')

## Base de Dados d4
# Bell d4
hist(eps_bell_d4, breaks = H, probability = TRUE, main = "PIT Bell (n = 1000)", col = cores[4])
abline(h = 1, lty = 2, col = 'red')
# Poisson d4
hist(eps_pois_d4, breaks = H, probability = TRUE, main = "PIT Poisson (n = 1000)", col = cores[4])
abline(h = 1, lty = 2, col = 'red')


##### Q-Q Plots ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Q-Q Plot Bell d1
qqplot(x = qunif(ppoints(length(eps_bell_d1))), y = eps_bell_d1, 
       main = "Q-Q Plot Bell (n = 50)", pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
# Q-Q Plot Poisson d1
qqplot(x = qunif(ppoints(length(eps_pois_d1))), y = eps_pois_d1, 
       main = "Q-Q Plot Poisson (n = 50)", pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

## Base de Dados d2
# Q-Q Plot Bell d2
qqplot(x = qunif(ppoints(length(eps_bell_d2))), y = eps_bell_d2, 
       main = "Q-Q Plot Bell (n = 100)", col = cores[2], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

# Q-Q Plot Poisson d2
qqplot(x = qunif(ppoints(length(eps_pois_d2))), y = eps_pois_d2, 
       main = "Q-Q Plot Poisson (n = 100)", col = cores[2], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

## Base de Dados d3
# Q-Q Plot Bell d3
qqplot(x = qunif(ppoints(length(eps_bell_d3))), y = eps_bell_d3,
       main = "Q-Q Plot Bell (n = 500)", col = cores[3], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

# Q-Q Plot Poisson d3
qqplot(x = qunif(ppoints(length(eps_pois_d3))), y = eps_pois_d3,
       main = "Q-Q Plot Poisson (n = 500)", col = cores[3], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

## Base de Dados d4
# Q-Q Plot Bell d4
qqplot(x = qunif(ppoints(length(eps_bell_d4))), y = eps_bell_d4,
       main = "Q-Q Plot Bell (n = 1000)", col = cores[4], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

# Q-Q Plot Poisson d4
qqplot(x = qunif(ppoints(length(eps_pois_d4))), y = eps_pois_d4,
       main = "Q-Q Plot Poisson (n = 1000)", col = cores[4], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')


#### Normalidade dos Resíduos ----

## Base de Dados d1
# Bell d1
eps.norm_bell_d1 = qnorm(eps_bell_d1)
eps.norm_bell_d1 = eps.norm_bell_d1[is.finite(eps.norm_bell_d1)] # Remove Inf e -Inf
p_val_bell_d1 <- shapiro.test(eps.norm_bell_d1)$p.value
p_format_bell_d1 <- formatC(p_val_bell_d1, format = "f", digits = 4)
# Poisson d1
eps.norm_pois_d1 = qnorm(eps_pois_d1)
eps.norm_pois_d1 = eps.norm_pois_d1[is.finite(eps.norm_pois_d1)]
p_val_pois_d1 <- shapiro.test(eps.norm_pois_d1)$p.value
p_format_pois_d1 <- formatC(p_val_pois_d1, format = "f", digits = 4)

## Base de Dados d2
# Bell d2
eps.norm_bell_d2 = qnorm(eps_bell_d2)
eps.norm_bell_d2 = eps.norm_bell_d2[is.finite(eps.norm_bell_d2)]
p_val_bell_d2 <- shapiro.test(eps.norm_bell_d2)$p.value
p_format_bell_d2 <- formatC(p_val_bell_d2, format = "f", digits = 4)
# Poisson d2
eps.norm_pois_d2 = qnorm(eps_pois_d2)
eps.norm_pois_d2 = eps.norm_pois_d2[is.finite(eps.norm_pois_d2)]
p_val_pois_d2 <- shapiro.test(eps.norm_pois_d2)$p.value
p_format_pois_d2 <- formatC(p_val_pois_d2, format = "f", digits = 4)

## Base de Dados d3
# Bell d3
eps.norm_bell_d3 = qnorm(eps_bell_d3)
eps.norm_bell_d3 = eps.norm_bell_d3[is.finite(eps.norm_bell_d3)]
p_val_bell_d3 <- shapiro.test(eps.norm_bell_d3)$p.value
p_format_bell_d3 <- formatC(p_val_bell_d3, format = "f", digits = 4)
# Poisson d3
eps.norm_pois_d3 = qnorm(eps_pois_d3)
eps.norm_pois_d3 = eps.norm_pois_d3[is.finite(eps.norm_pois_d3)]
p_val_pois_d3 <- shapiro.test(eps.norm_pois_d3)$p.value
p_format_pois_d3 <- formatC(p_val_pois_d3, format = "f", digits = 4)

## Base de Dados d4
# Bell d4
eps.norm_bell_d4 = qnorm(eps_bell_d4)
eps.norm_bell_d4 = eps.norm_bell_d4[is.finite(eps.norm_bell_d4)]
p_val_bell_d4 <- shapiro.test(eps.norm_bell_d4)$p.value
p_format_bell_d4 <- formatC(p_val_bell_d4, format = "f", digits = 4)
# Poisson d4
eps.norm_pois_d4 = qnorm(eps_pois_d4)
eps.norm_pois_d4 = eps.norm_pois_d4[is.finite(eps.norm_pois_d4)]
p_val_pois_d4 <- shapiro.test(eps.norm_pois_d4)$p.value
p_format_pois_d4 <- formatC(p_val_pois_d4, format = "f", digits = 4)


##### Histogramas ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
hist(eps.norm_bell_d1, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 50)\n p-value =", p_format_bell_d1), 
     col = cores[1], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d1), sd = sd(eps.norm_bell_d1)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d1
hist(eps.norm_pois_d1, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 50)\n p-value =", p_format_pois_d1), 
     col = cores[1], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d1), sd = sd(eps.norm_pois_d1)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d2
# Bell d2
hist(eps.norm_bell_d2, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 100)\n p-value =", p_format_bell_d2), 
     col = cores[2], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d2), sd = sd(eps.norm_bell_d2)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d2
hist(eps.norm_pois_d2, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 100)\n p-value =", p_format_pois_d2), 
     col = cores[2], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d2), sd = sd(eps.norm_pois_d2)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d3
# Bell d3
hist(eps.norm_bell_d3, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 500)\n p-value =", p_format_bell_d3), 
     col = cores[3], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d3), sd = sd(eps.norm_bell_d3)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d3
hist(eps.norm_pois_d3, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 500)\n p-value =", p_format_pois_d3), 
     col = cores[3], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d3), sd = sd(eps.norm_pois_d3)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d4
# Bell d4
hist(eps.norm_bell_d4, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 1000)\n p-value =", p_format_bell_d4), 
     col = cores[4], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d4), sd = sd(eps.norm_bell_d4)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d4
hist(eps.norm_pois_d4, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 1000)\n p-value =", p_format_pois_d4), 
     col = cores[4], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d4), sd = sd(eps.norm_pois_d4)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)


##### Q-Q Plots ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
qqnorm(eps.norm_bell_d1, 
       main = "Q-Q Plot Normal\n Bell (n = 50)",
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d1, col = "red")
# Poisson d1
qqnorm(eps.norm_pois_d1, 
       main = "Q-Q Plot Normal\n Poisson (n = 50)",
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d1, col = "red")

## Base de Dados d2
# Bell d2
qqnorm(eps.norm_bell_d2, 
       main = "Q-Q Plot Normal\n Bell (n = 100)", col = cores[2],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d2, col = "red")
# Poisson d2
qqnorm(eps.norm_pois_d2, 
       main = "Q-Q Plot Normal\n Poisson (n = 100)", col = cores[2],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d2, col = "red")

## Base de Dados d3
# Bell d3
qqnorm(eps.norm_bell_d3, 
       main = "Q-Q Plot Normal\n Bell (n = 500)", col = cores[3],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d3, col = "red")
# Poisson d3
qqnorm(eps.norm_pois_d3, 
       main = "Q-Q Plot Normal\n Poisson (n = 500)", col = cores[3],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d3, col = "red")

## Base de Dados d4
# Bell d4
qqnorm(eps.norm_bell_d4, 
       main = "Q-Q Plot Normal\n Bell (n = 1000)", col = cores[4],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d4, col = "red")
#Poisson d4
qqnorm(eps.norm_pois_d4, 
       main = "Q-Q Plot Normal\n Poisson (n = 1000)", col = cores[4],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d4, col = "red")


rm(list = setdiff(ls(), c("scaled_residuals_bell", "scaled_residuals_pois",
                          "scaled_residuals_nb", "H")))
## Variável Explicativa Exponencial ----
# Função de dados com variável explicativa exponencial
dadosunif <- function(n) {
  mX = matrix(NA, nrow = n, ncol = 2)
  vBeta = matrix(c(1, 2))
  vY = rep(NA, n)
  vMu = rep(NA, n)
  vX1 = rep(1, n)
  vX2_raw = runif(n, min = 0, max = 1)
  vX2_scaled = scale(vX2_raw, 
                     scale = FALSE, center = TRUE) # padronizando vX2
  for(i in 1:n){
    mX[i, 1] = vX1[i]
    mX[i, 2] = vX2_scaled[i]
    vMu[i] = exp(t(mX[i, ]) %*% vBeta)
    vY[i] = rbell(n = 1, theta = LambertW::W(vMu[i]))
  }
  return(list(vY = vY, mX = mX ))
}

dadosexp <- function(n) {
  mX = matrix(NA, nrow = n, ncol = 2)
  vBeta = matrix(c(1, 2))
  vY = rep(NA, n)
  vMu = rep(NA, n)
  vX1 = rep(1, n)
  vX2_raw = rexp(n, rate = 5)
  vX2_scaled = scale(vX2_raw, 
                     scale = FALSE, center = TRUE) # padronizando vX2
  for(i in 1:n){
    mX[i, 1] = vX1[i]
    mX[i, 2] = vX2_scaled[i]
    vMu[i] = exp(t(mX[i, ]) %*% vBeta)
    vY[i] = rbell(n = 1, theta = LambertW::W(vMu[i]))
  }
  return(list(vY = vY, mX = mX ))
}
# Vetor de cores usados
cores <- c("#BDD7E7", "#BAE4B3", "#2171B5", "#238B45")

set.seed(20)
d1 <- dadosunif(100)
d2 <- dadosexp(100)
d3 <- dadosunif(1000)
d4 <- dadosexp(1000)


### Visualização ----
par(mfrow = c(2, 2))
datasets <- list(d1, d2, d3, d4)
n_obs <- c(100, 100, 1000, 1000)
d_obs <- c("Uniforme", "Exponencial", "Uniforme", "Exponencial")

for (i in 1:4) {
  hist(datasets[[i]]$vY,
       prob = TRUE,
       main = paste("Histograma de vY\n Explicativa", d_obs[i], "\n(n =", n_obs[i], ")"),
       col = cores[i],
       xlab = "vY",
       ylab = "Densidade")
  
  lines(density(datasets[[i]]$vY))
}

### Modelos MLG ----
## Base de Dados d1
# Bell d1
out.belld1 = glm(vY ~ mX[, 2], data = d1,
                 family = bell(link = log))
# Poisson d1
out.poisd1 = glm(vY ~ mX[, 2], data = d1,
                 family = poisson(link = log))
## Base de Dados d2
# Bell d2
out.belld2 = glm(vY ~ mX[, 2], data = d2,
                 family = bell(link = log))
# Poisson d2
out.poisd2 = glm(vY ~ mX[, 2], data = d2,
                 family = poisson(link = log))
## Base de Dados d3
# Bell d3
out.belld3 = glm(vY ~ mX[, 2], data = d3,
                 family = bell(link = log))
# Poisson d3
out.poisd3 = glm(vY ~ mX[, 2], data = d3,
                 family = poisson(link = log))
## Base de Dados d4
# Bell d4
out.belld4 = glm(vY ~ mX[, 2], data = d4,
                 family = bell(link = log))
# Poisson d4
out.poisd4 = glm(vY ~ mX[, 2], data = d4,
                 family = poisson(link = log))


### AIC ----
# Comparando AIC dos modelos
tabela_aic <- data.frame(
  Amostra = c("Unif n = 100", "Exp n = 100", "Unif n = 1000", "Exp n = 1000"),
  Bell    = c(round(AIC(out.belld1), 2),
              round(AIC(out.belld2), 2),
              round(AIC(out.belld3), 2),
              round(AIC(out.belld4), 2)),
  
  Poisson = c(round(AIC(out.poisd1), 2),
              round(AIC(out.poisd2), 2),
              round(AIC(out.poisd3), 2),
              round(AIC(out.poisd4), 2))
)

kable(tabela_aic,
      col.names = c("Amostra (n)", "Modelo Bell", "Modelo Poisson"),
      align = c("l", "c", "c"),
      row.names = FALSE,
      caption = "Comparação do AIC dos Modelos Bell e Poisson por Amostra")


### Estimativa de Betas ----
# n = 100 Uniforme
summary_belld1 <- summary(out.belld1)
coef_belld1 <- summary_belld1$coefficients[, "Estimate"]

summary_poisd1 <- summary(out.poisd1)
coef_poisd1 <- summary_poisd1$coefficients[, "Estimate"]

# n = 100 Exponencial
summary_belld2 <- summary(out.belld2)
coef_belld2 <- summary_belld2$coefficients[, "Estimate"]

summary_poisd2 <- summary(out.poisd2)
coef_poisd2 <- summary_poisd2$coefficients[, "Estimate"]

# n = 1000 Uniforme
summary_belld3 <- summary(out.belld3)
coef_belld3 <- summary_belld3$coefficients[, "Estimate"]

summary_poisd3 <- summary(out.poisd3)
coef_poisd3 <- summary_poisd3$coefficients[, "Estimate"]

# n = 1000 Exponencial
summary_belld4 <- summary(out.belld4)
coef_belld4 <- summary_belld4$coefficients[, "Estimate"]

summary_poisd4 <- summary(out.poisd4)
coef_poisd4 <- summary_poisd4$coefficients[, "Estimate"]

tabela_betas <- data.frame(
  Amostra = c("Unif n = 100", "Exp n = 100", "Unif n = 1000", "Exp n = 1000"),
  B0_Bell = c(round(coef_belld1[1], 2),
              round(coef_belld2[1], 2),
              round(coef_belld3[1], 2),
              round(coef_belld4[1], 2)),
  
  B0_Pois = c(round(coef_poisd1[1], 2),
              round(coef_poisd2[1], 2),
              round(coef_poisd3[1], 2),
              round(coef_poisd4[1], 2)),
  
  B1_Bell = c(round(coef_belld1[2], 2),
              round(coef_belld2[2], 2),
              round(coef_belld3[2], 2),
              round(coef_belld4[2], 2)),
  
  B1_Pois = c(round(coef_poisd1[2], 2),
              round(coef_poisd2[2], 2),
              round(coef_poisd3[2], 2),
              round(coef_poisd4[2], 2))
)

# Gerando a tabela de Betas
kable(tabela_betas,
      col.names = c("Amostra (n)", 
                    "Beta 0 (Bell)", "Beta 0 (Poisson)",
                    "Beta 1 (Bell)", "Beta 1 (Poisson)"),
      align = c("l", "c", "c", "c", "c"),
      row.names = FALSE,
      caption = "Estimativas dos Parâmetros Beta 0 e Beta 1 por Amostra")



### Análise de Resíduos ----

#### Resíduos PIT ----
set.seed(20)
## Base de Dados d1
# Bell d1
L_bell_d1 = scaled_residuals_bell(out.belld1, random = TRUE, B = 250)
eps_bell_d1 = L_bell_d1$scaled_residuals
# Poisson d1
L_pois_d1 = scaled_residuals_pois(out.poisd1, random = TRUE, B = 250)
eps_pois_d1 = L_pois_d1$scaled_residuals

## Base de Dados d2
# Bell d2
L_bell_d2 = scaled_residuals_bell(out.belld2, random = TRUE, B = 250)
eps_bell_d2 = L_bell_d2$scaled_residuals
# Poisson d2
L_pois_d2 = scaled_residuals_pois(out.poisd2, random = TRUE, B = 250)
eps_pois_d2 = L_pois_d2$scaled_residuals

## Base de Dados d3
# Bell d3
L_bell_d3 = scaled_residuals_bell(out.belld3, random = TRUE, B = 250)
eps_bell_d3 = L_bell_d3$scaled_residuals
# Poisson d3
L_pois_d3 = scaled_residuals_pois(out.poisd3, random = TRUE, B = 250)
eps_pois_d3 = L_pois_d3$scaled_residuals

## Base de Dados d4
# Bell d4
L_bell_d4 = scaled_residuals_bell(out.belld4, random = TRUE, B = 250)
eps_bell_d4 = L_bell_d4$scaled_residuals
# Poisson d4
L_pois_d4 = scaled_residuals_pois(out.poisd4, random = TRUE, B = 250)
eps_pois_d4 = L_pois_d4$scaled_residuals


##### Histogramas ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
hist(eps_bell_d1, breaks = H, probability = TRUE, main = "PIT Bell (n = 100)", col = cores[1])
abline(h = 1, lty = 2, col = 'red')
# Poisson d1
hist(eps_pois_d1, breaks = H, probability = TRUE, main = "PIT Poisson (n = 100)", col = cores[1])
abline(h = 1, lty = 2, col = 'red')
## Base de Dados d2
# Bell d2
hist(eps_bell_d2, breaks = H, probability = TRUE, main = "PIT Bell (n = 100)", col = cores[2])
abline(h = 1, lty = 2, col = 'red')
# Poisson d2
hist(eps_pois_d2, breaks = H, probability = TRUE, main = "PIT Poisson (n = 100)", col = cores[2])
abline(h = 1, lty = 2, col = 'red')
## Base de Dados d3
# Bell d3
hist(eps_bell_d3, breaks = H, probability = TRUE, main = "PIT Bell (n = 1000)", col = cores[3])
abline(h = 1, lty = 2, col = 'red')
# Poisson d3
hist(eps_pois_d3, breaks = H, probability = TRUE, main = "PIT Poisson (n = 1000)", col = cores[3])
abline(h = 1, lty = 2, col = 'red')
## Base de Dados d4
# Bell d4
hist(eps_bell_d4, breaks = H, probability = TRUE, main = "PIT Bell (n = 1000)", col = cores[4])
abline(h = 1, lty = 2, col = 'red')
# Poisson d4
hist(eps_pois_d4, breaks = H, probability = TRUE, main = "PIT Poisson (n = 1000)", col = cores[4])
abline(h = 1, lty = 2, col = 'red')


##### Q-Q Plots ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
qqplot(x = qunif(ppoints(length(eps_bell_d1))), y = eps_bell_d1, 
       main = "Q-Q Plot Bell (n = 100)", col = cores[1], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
# Poisson d1
qqplot(x = qunif(ppoints(length(eps_pois_d1))), y = eps_pois_d1, 
       main = "Q-Q Plot Poisson (n = 100)", col = cores[1], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
## Base de Dados d2
# Bell d2
qqplot(x = qunif(ppoints(length(eps_bell_d2))), y = eps_bell_d2, 
       main = "Q-Q Plot Bell (n = 100)", col = cores[2], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
# Poisson d2
qqplot(x = qunif(ppoints(length(eps_pois_d2))), y = eps_pois_d2, 
       main = "Q-Q Plot Poisson (n = 100)", col = cores[2], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
## Base de Dados d3
# Bell d3
qqplot(x = qunif(ppoints(length(eps_bell_d3))), y = eps_bell_d3, 
       main = "Q-Q Plot Bell (n = 1000)", col = cores[3], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
# Poisson d3
qqplot(x = qunif(ppoints(length(eps_pois_d3))), y = eps_pois_d3, 
       main = "Q-Q Plot Poisson (n = 1000)", col = cores[3], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
## Base de Dados d4
# Bell d4
qqplot(x = qunif(ppoints(length(eps_bell_d4))), y = eps_bell_d4, 
       main = "Q-Q Plot Bell (n = 1000)", col = cores[4], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')
# Poisson d4
qqplot(x = qunif(ppoints(length(eps_pois_d4))), y = eps_pois_d4, 
       main = "Q-Q Plot Poisson (n = 1000)", col = cores[4], pch = 16,
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1, lty = 2, col = 'red')

#### Normalidade dos Resíduos ----

## Base de Dados d1
# Bell d1
eps.norm_bell_d1 = qnorm(eps_bell_d1)
eps.norm_bell_d1 = eps.norm_bell_d1[is.finite(eps.norm_bell_d1)]
p_val_bell_d1 <- shapiro.test(eps.norm_bell_d1)$p.value
p_format_bell_d1 <- formatC(p_val_bell_d1, format = "f", digits = 4)
# Poisson d1
eps.norm_pois_d1 = qnorm(eps_pois_d1)
eps.norm_pois_d1 = eps.norm_pois_d1[is.finite(eps.norm_pois_d1)]
p_val_pois_d1 <- shapiro.test(eps.norm_pois_d1)$p.value
p_format_pois_d1 <- formatC(p_val_pois_d1, format = "f", digits = 4)

## Base de Dados d2
# Bell d2
eps.norm_bell_d2 = qnorm(eps_bell_d2)
eps.norm_bell_d2 = eps.norm_bell_d2[is.finite(eps.norm_bell_d2)]
p_val_bell_d2 <- shapiro.test(eps.norm_bell_d2)$p.value
p_format_bell_d2 <- formatC(p_val_bell_d2, format = "f", digits = 4)
# Poisson d2
eps.norm_pois_d2 = qnorm(eps_pois_d2)
eps.norm_pois_d2 = eps.norm_pois_d2[is.finite(eps.norm_pois_d2)]
p_val_pois_d2 <- shapiro.test(eps.norm_pois_d2)$p.value
p_format_pois_d2 <- formatC(p_val_pois_d2, format = "f", digits = 4)

## Base de Dados d3
# Bell d3
eps.norm_bell_d3 = qnorm(eps_bell_d3)
eps.norm_bell_d3 = eps.norm_bell_d3[is.finite(eps.norm_bell_d3)]
p_val_bell_d3 <- shapiro.test(eps.norm_bell_d3)$p.value
p_format_bell_d3 <- formatC(p_val_bell_d3, format = "f", digits = 4)
# Poisson d3
eps.norm_pois_d3 = qnorm(eps_pois_d3)
eps.norm_pois_d3 = eps.norm_pois_d3[is.finite(eps.norm_pois_d3)]
p_val_pois_d3 <- shapiro.test(eps.norm_pois_d3)$p.value
p_format_pois_d3 <- formatC(p_val_pois_d3, format = "f", digits = 4)

## Base de Dados d4
# Bell d4
eps.norm_bell_d4 = qnorm(eps_bell_d4)
eps.norm_bell_d4 = eps.norm_bell_d4[is.finite(eps.norm_bell_d4)]
p_val_bell_d4 <- shapiro.test(eps.norm_bell_d4)$p.value
p_format_bell_d4 <- formatC(p_val_bell_d4, format = "f", digits = 4)
# Poisson d4
eps.norm_pois_d4 = qnorm(eps_pois_d4)
eps.norm_pois_d4 = eps.norm_pois_d4[is.finite(eps.norm_pois_d4)]
p_val_pois_d4 <- shapiro.test(eps.norm_pois_d4)$p.value
p_format_pois_d4 <- formatC(p_val_pois_d4, format = "f", digits = 4)


##### Histogramas ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
hist(eps.norm_bell_d1, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 100)\n p-value =", p_format_bell_d1), 
     col = cores[1], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d1), sd = sd(eps.norm_bell_d1)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d1
hist(eps.norm_pois_d1, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 100)\n p-value =", p_format_pois_d1), 
     col = cores[1], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d1), sd = sd(eps.norm_pois_d1)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d2
# Bell d2
hist(eps.norm_bell_d2, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 100)\n p-value =", p_format_bell_d2), 
     col = cores[2], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d2), sd = sd(eps.norm_bell_d2)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d2
hist(eps.norm_pois_d2, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 100)\n p-value =", p_format_pois_d2), 
     col = cores[2], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d2), sd = sd(eps.norm_pois_d2)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d3
# Bell d3
hist(eps.norm_bell_d3, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 1000)\n p-value =", p_format_bell_d3), 
     col = cores[3], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d3), sd = sd(eps.norm_bell_d3)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d3
hist(eps.norm_pois_d3, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 1000)\n p-value =", p_format_pois_d3), 
     col = cores[3], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d3), sd = sd(eps.norm_pois_d3)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)

## Base de Dados d4
# Bell d4
hist(eps.norm_bell_d4, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Bell (n = 1000)\n p-value =", p_format_bell_d4), 
     col = cores[4], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_bell_d4), sd = sd(eps.norm_bell_d4)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)
# Poisson d4
hist(eps.norm_pois_d4, probability = TRUE, 
     main = paste("Resíduos Normalizados\n Poisson (n = 1000)\n p-value =", p_format_pois_d4), 
     col = cores[4], xlab = "Resíduos", ylab = "Densidade")
curve(dnorm(x, mean = mean(eps.norm_pois_d4), sd = sd(eps.norm_pois_d4)),  
      add = TRUE, col = "red", lty = 2, lwd = 2)


##### Q-Q Plots ----
par(mfrow = c(1, 2))
## Base de Dados d1
# Bell d1
qqnorm(eps.norm_bell_d1, 
       main = "Q-Q Plot Normal\n Bell (n = 100)", col = cores[1],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d1, col = "red")
# Poisson d1
qqnorm(eps.norm_pois_d1, 
       main = "Q-Q Plot Normal\n Poisson (n = 100)", col = cores[1],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d1, col = "red")

## Base de Dados d2
# Bell d2
qqnorm(eps.norm_bell_d2, 
       main = "Q-Q Plot Normal\n Bell (n = 100)", col = cores[2],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d2, col = "red")
# Poisson d2
qqnorm(eps.norm_pois_d2, 
       main = "Q-Q Plot Normal\n Poisson (n = 100)", col = cores[2],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d2, col = "red")

## Base de Dados d3
# Bell d3
qqnorm(eps.norm_bell_d3, 
       main = "Q-Q Plot Normal\n Bell (n = 1000)", col = cores[3],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d3, col = "red")
# Poisson d3
qqnorm(eps.norm_pois_d3, 
       main = "Q-Q Plot Normal\n Poisson (n = 1000)", col = cores[3],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d3, col = "red")

## Base de Dados d4
# Bell d4
qqnorm(eps.norm_bell_d4, 
       main = "Q-Q Plot Normal\n Bell (n = 1000)", col = cores[4],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_bell_d4, col = "red")
# Poisson d4
qqnorm(eps.norm_pois_d4, 
       main = "Q-Q Plot Normal\n Poisson (n = 1000)", col = cores[4],
       xlab = "Quantis Teóricos", ylab = "Quantis Amostrais")
qqline(eps.norm_pois_d4, col = "red")




rm(list = setdiff(ls(), c("scaled_residuals_bell", "scaled_residuals_pois",
                          "scaled_residuals_nb", "H")))
# Dados Reais ----

## Base de Dados ----
# https://www.kaggle.com/datasets/adilshamim8/math-students
data <- read.csv("~/Math-Students.csv")
## Definindo variáveis como fatores para modelo final
# Variáveis dos Modelos Bell e Binomial Negativa
data$studytime <- as.factor(data$studytime)
data$health <- as.factor(data$health)
# Variáveis do Modelo Poisson
data$famsize <- as.factor(data$famsize)
data$Medu <- as.factor(data$Medu)
data$traveltime <- as.factor(data$traveltime)
data$freetime <- as.factor(data$freetime)
data$goout <- as.factor(data$goout)
data$Walc <- as.factor(data$Walc)
# Nota Total
data$GT <- data$G1 + data$G2 + data$G3 #Grade Total
# Tamanho Amostra
n = nrow(data)

## Análise Descritiva ----
### Sumário de Nota Total ----
# Variável G1
media_g1    <- mean(data$G1)
mediana_g1  <- median(data$G1)
var_g1      <- var(data$G1)
q1_g1       <- quantile(data$G1, 0.25)
q3_g1       <- quantile(data$G1, 0.75)
min_g1      <- min(data$G1)
max_g1      <- max(data$G1)
# Variável G2
media_g2    <- mean(data$G2)
mediana_g2  <- median(data$G2)
var_g2      <- var(data$G2)
q1_g2       <- quantile(data$G2, 0.25)
q3_g2       <- quantile(data$G2, 0.75)
min_g2      <- min(data$G2)
max_g2      <- max(data$G2)
# Variável G3
media_g3    <- mean(data$G3)
mediana_g3  <- median(data$G3)
var_g3      <- var(data$G3)
q1_g3       <- quantile(data$G3, 0.25)
q3_g3       <- quantile(data$G3, 0.75)
min_g3      <- min(data$G3)
max_g3      <- max(data$G3)
# Variável G4 (GT)
media_gt    <- mean(data$GT)
mediana_gt  <- median(data$GT)
var_gt      <- var(data$GT)
q1_gt       <- quantile(data$GT, 0.25)
q3_gt       <- quantile(data$GT, 0.75)
min_gt      <- min(data$GT)
max_gt      <- max(data$GT)

tabela_sumario <- data.frame(
  Variavel = c("Nota 1", "Nota 2", "Nota 3", "Nota Total"),
  
  Média = c(round(media_g1, 2), 
            round(media_g2, 2), 
            round(media_g3, 2), 
            round(media_gt, 2)),
  
  Variancia = c(round(var_g1, 2), 
                round(var_g2, 2), 
                round(var_g3, 2), 
                round(var_gt, 2)),
  
  Minimo = c(round(min_g1, 2),
             round(min_g2, 2),
             round(min_g3, 2),
             round(min_gt, 2)),
  
  Q1 = c(round(q1_g1, 2), 
         round(q1_g2, 2), 
         round(q1_g3, 2), 
         round(q1_gt, 2)),
  
  Mediana = c(round(mediana_g1, 2), 
              round(mediana_g2, 2), 
              round(mediana_g3, 2), 
              round(mediana_gt, 2)),
  
  Q3 = c(round(q3_g1, 2), 
         round(q3_g2, 2), 
         round(q3_g3, 2), 
         round(q3_gt, 2)),
  
  Maximo = c(round(max_g1, 2),
             round(max_g2, 2),
             round(max_g3, 2),
             round(max_gt, 2))
)

kable(tabela_sumario,
      col.names = c("Variável", "Média", "Var", "Min", "Q1", "Mediana", "Q3", "Max"),
      align = c("l", "c", "c", "c", "c", "c", "c", "c"),
      row.names = FALSE,
      caption = "Sumário das Notas 1, 2, 3 e Total")

### Densidade de Nota Total ----
par(mfrow = c(1, 1))
ggplot(data) +
  geom_density(aes(x = G1, color = "G1"), linewidth = 1.2) +
  geom_density(aes(x = G2, color = "G2"), linewidth = 1.2) +
  geom_density(aes(x = G3, color = "G3"), linewidth = 1.2) +
  geom_density(aes(x = GT, color = "GT"), linewidth = 1.2) +
  
  # Aplica a paleta do RColorBrewer manualmente para as 4 linhas
  scale_color_manual(values = brewer.pal(4, "Set1")) + 
  
  labs(
    title = "Comparação de Densidade das Notas",
    x = "Valores",
    y = "Densidade",
    color = "Provas"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
  )

### Sexo ----
ggplot(data, aes(x = sex, y = GT, fill = sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink", "lightblue"), guide = "none") +
  scale_x_discrete(labels = c("F" = "Feminino", "M" = "Masculino")) +
  labs(
    title = "Distribuição das Notas por Sexo",
    x = "Sexo",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()

### Tempo de Estudo Semanal ----
ggplot(data, aes(x = studytime, y = GT, fill = studytime)) +
  geom_boxplot() +
  scale_fill_manual(values = brewer.pal(4, "Purples"), guide = "none") +
  scale_x_discrete(
    labels =c("1" = "Até 15 minutos", "2" = "15 a 30 minutos",
              "3" = "30 minutos a 1 hora", "4" = "mais de uma hora")) +
  labs(
    title = "Distribuição das Notas por Tempo de Estudo",
    x = "Tempo de Estudo",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()

### Repetências ----
ggplot(data, aes(x = as.factor(failures), y = GT, fill = as.factor(failures))) +
  geom_boxplot() +
  scale_fill_manual(values = brewer.pal(4, "Reds"), guide = "none") +
  labs(
    title = "Distribuição das Notas por Repetências",
    x = "Repetências",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()

### Suporte Extra Escolar ----
ggplot(data, aes(x = schoolsup, y = GT, fill = schoolsup)) +
  geom_boxplot() +
  scale_fill_manual(values = c("indianred1", "lightblue1"), guide = "none") +
  scale_x_discrete(labels = c("no" = "Não", "yes" = "Sim")) +
  labs(
    title = "Distribuição das Notas por Suporte Extra Escolar",
    x = "Recebe Suporte Educacional Adicional da Escola",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()

### Ensino Superior ----
ggplot(data, aes(x = higher, y = GT, fill = higher)) +
  geom_boxplot() +
  scale_fill_manual(values = c("indianred1", "lightblue1"), guide = "none") +
  scale_x_discrete(labels = c("no" = "Não", "yes" = "Sim")) +
  labs(
    title = "Distribuição das Notas Entre Aqueles que \nAlmejam ou Não Ensino Superior",
    x = "Almeja Ensino Superior",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()

### Saúde ----
ggplot(data, aes(x = health, y = GT, fill = health)) +
  geom_boxplot() +
  scale_fill_manual(values = brewer.pal(5, "RdYlGn"), guide = "none") +
  scale_x_discrete(
    labels = c("1" = "Muito Ruim", "2" = "Ruim", "3" = "Mediano", "4" = "Bom", "5" = "Muito Bom")) +
  labs(
    title = "Distribuição das Notas pelo Estado de Saúde",
    x = "Estado de Saúde",
    y = "Nota Total (GT)"
  ) +
  theme_minimal()


## Modelos ----
cores <- brewer.pal(3, "Set2")
out.bell = glm(GT ~ sex + studytime + failures + schoolsup + higher + health,
               family = bell(link = log), data = data)

out.pois = glm(GT ~ sex + famsize + Medu + Mjob + Fjob + traveltime + studytime +
                 failures + schoolsup + famsup + higher + internet + romantic +
                 freetime + goout + Walc + health + absences,
               family = poisson(link = log), data = data)

out.nbin = glm.nb(GT ~ sex + studytime + failures + schoolsup +
                    higher + health, link = log, data = data)


### AIC ----
tf <- data.frame(
  Metrica = c(
    "AIC"
  ),
  Modelo_Bell = c(
    round(AIC(out.bell), 2)
  ),
  Modelo_Poisson = c(
    round(AIC(out.pois), 2)
  ),
  Modelo_BinomialN = c(
    round(AIC(out.nbin), 2)
  )
)
kable(tf,
      col.names = c("Métrica / Parâmetro", "Bell", "Poisson", "Binomial Negativa"),
      align = c("l", "c", "c", "c"),
      caption = "Comparação do AIC dos Modelos")


### Resíduos PIT ----
set.seed(20)

L_bell = scaled_residuals_bell(out.bell, random = TRUE, B = 250)
eps_bell = L_bell$scaled_residuals

L_pois = scaled_residuals_pois(out.pois, random = TRUE, B = 250)
eps_pois = L_pois$scaled_residuals

L_nbin = scaled_residuals_nb(out.nbin, random = TRUE, B = 250)
eps_nbin = L_nbin$scaled_residuals


#### Histogramas ----
# Bell
hist(eps_bell, breaks = H, probability = TRUE,
     col = cores[1], main = "PIT Bell")
abline(h = 1, lty = 2, col = 'red') # referência
# Poisson
hist(eps_pois, breaks = H, probability = TRUE,
     col = cores[2], main = "PIT Poisson")
abline(h = 1, lty = 2, col = 'red')
# Binomial Negativa
hist(eps_nbin, breaks = H, probability = TRUE,
     col = cores[3], main = "PIT Binomial Negativa")
abline(h = 1, lty = 2, col = 'red')


#### Q-Q Plots ----
# Bell
qqplot(x = qunif(ppoints(n)), y = eps_bell,
       col = cores[1], main = "Q-Q Plot Bell",
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1)
# Poisson
qqplot(x = qunif(ppoints(n)), y = eps_pois,
       col = cores[2], main = "Q-Q Plot Poisson",
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1)
# Binomial Negativa
qqplot(x = qunif(ppoints(n)), y = eps_nbin,
       col = cores[3], main = "Q-Q Plot Binomial Negativa",
       xlab = "Quantis Teóricos Uniforme", ylab = "Quantis Amostrais")
abline(a = 0, b = 1)


### Normalidade ----
# Bell
eps.norm_bell = qnorm(eps_bell)
eps.norm_bell = eps.norm_bell[is.finite(eps.norm_bell)]
pbell <- shapiro.test(eps.norm_bell)$p.value
format_pbell <- formatC(pbell, format = "f", digits = 4)
# Poisson
eps.norm_pois = qnorm(eps_pois)
eps.norm_pois = eps.norm_pois[is.finite(eps.norm_pois)]
ppois <- shapiro.test(eps.norm_pois)$p.value
format_ppois <- formatC(ppois, format = "f", digits = 4)
# Binomial Negativa
eps.norm_nbin = qnorm(eps_nbin)
eps.norm_nbin = eps.norm_nbin[is.finite(eps.norm_nbin)]
pnbin <- shapiro.test(eps.norm_nbin)$p.value
format_pnbin <- formatC(pnbin, format = "f", digits = 4)

#### Histogramas ----
# Bell
hist(eps.norm_bell, probability = TRUE, col = cores[1], main =
       paste("Resíduos Normalizados Bell \n p-value =", format_pbell))
curve(dnorm(x, mean = mean(eps.norm_bell), sd = sd(eps.norm_bell)),
      add = TRUE, col = "red", lty = 2)
# Poisson
hist(eps.norm_pois, probability = TRUE, col = cores[2], main =
       paste("Resíduos Normalizados Poisson \n p-value =", format_ppois))
curve(dnorm(x, mean = mean(eps.norm_pois), sd = sd(eps.norm_pois)),
      add = TRUE, col = "red", lty = 2)
# Binomial Negativa
hist(eps.norm_nbin, probability = TRUE, col = cores[3], main =
       paste("Resíduos Normalizados Binomial Negativa \n p-value =", format_pnbin))
curve(dnorm(x, mean = mean(eps.norm_nbin), sd = sd(eps.norm_nbin)),
      add = TRUE, col = "red", lty = 2)

#### Q-Q Plots ----
# Bell
qqnorm(eps.norm_bell, 
       col = cores[1], main = "Q-Q Plot Bell",
       xlab = "Quantis Teóricos Normal", ylab = "Quantis Amostrais")
qqline(eps.norm_bell)
# Poisson
qqnorm(eps.norm_pois, 
       col = cores[2], main = "Q-Q Plot Poisson",
       xlab = "Quantis Teóricos Normal", ylab = "Quantis Amostrais")
qqline(eps.norm_pois)
# Binomial Negativa
qqnorm(eps.norm_nbin, 
       col = cores[3], main = "Q-Q Plot Binomial Negativa",
       xlab = "Quantis Teóricos Normal", ylab = "Quantis Amostrais")
qqline(eps.norm_nbin)

## Modelo Final ----
coef_summary <- summary(out.bell)$coefficients
estimativas <- coef_summary[, 1]
p_valores <- coef_summary[, 4]

# Intervalo de Confiança (IC 95%) do modelo
ic_modelo <- confint(out.bell)

# Razão de Chances e seu IC
razao_chances <- exp(estimativas)
ic_exponencializado <- exp(ic_modelo)

tabela_rc <- data.frame(
  Razao_Chances = razao_chances,
  IC_2.5 = ic_exponencializado[, 1],
  IC_97.5 = ic_exponencializado[, 2],
  p_valor = p_valores
)

kable(
  tabela_rc,
  digits = 4,
  col.names = c("Risco Relativo", "2.5 %", "97.5 %", "p-valor"),
  align = c("c", "c", "c", "c"),
  caption = "Resultados do MLG Bell"
)
