# FULL MINIMAL EXAMPLE (outside unmarked) --------------------------------------

library(future)
library(future.apply)

plan(multisession, workers = parallel::detectCores())

nb_rep <- 100
res_fit <- future_lapply(seq_len(nb_rep), function(id_rep) {
  ## Simulation parameters
  S <- 200
  J_i <- sample(x = 1, size = S, replace = TRUE)
  alpha0 <- 0.5
  alpha1 <- 1
  beta0 <- -1
  beta1 <- 1

  ## Site occupancy covariate
  siteCov <- rnorm(S)

  ## Occupancy
  psi <- plogis(alpha0 + alpha1 * siteCov)
  Z <- rbinom(S, 1, psi)

  ## Detection covariate
  obsCovFunctions <- vector("list", length = S)
  for (i in 1:S) {
    obsCovFunctions[[i]] <- vector("list", length = J_i[i])
    for (j in 1:J_i[i]) {
      eval(parse(text = paste(
        "obsCovFunctions[[i]][[j]] <- function(t) {",
        rnorm(1, sd = 0.1), "+", rnorm(1, sd = 0.1), "* t }"
      )))
    }
  }

  ## Simulation of detection events
  ovspg_lambda <- max(sapply(seq_along(obsCovFunctions), \(x) {
    max(sapply(seq_along(obsCovFunctions[[x]]), \(y) {
      max(sapply(seq(0, 1, length.out = 500), \(t) {
        exp(beta0 + beta1 * obsCovFunctions[[x]][[y]](t))
      }))
    }))
  })) * 1.05

  detec_list <- vector("list", S)
  T0 <- 0
  T1 <- 1
  for (i in 1:S) {
    detec_list[[i]] <- vector("list", J_i[i])
    for (j in 1:J_i[i]) {
      if (Z[i] == 1) {
        ovspg_times <- sort(runif(
          n = rpois(n = 1, lambda = ovspg_lambda * (T1 - T0)),
          min = T0, max = T1
        ))
        ovspg_exp_lambda <- exp(beta0 + beta1 * obsCovFunctions[[i]][[j]](ovspg_times))
        detec_list[[i]][[j]] <- ovspg_times[
          as.logical(rbinom(length(ovspg_times), 1, ovspg_exp_lambda / ovspg_lambda))
        ]
      } else {
        detec_list[[i]][[j]] <- numeric(0)
      }
    }
  }
  ndet <- sum(sapply(detec_list, \(x){
    sum(lengths(x))
  }))

  ## Plot
  if (F) {
    library(ggplot2)
    print(
      ggplot() +
        facet_grid(j ~ i) +
        geom_rect(
          data = data.frame(
            i = unlist(sapply(seq_len(S), \(i){
              rep(i, J_i[i])
            })),
            j = unlist(sapply(J_i, seq_len)),
            Z = unlist(sapply(seq_len(S), \(i){
              rep(as.logical(Z[i]), J_i[i])
            }))
          ),
          aes(
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = Z
          )
        ) +
        geom_line(
          data = bind_rows(
            lapply(1:S, \(i){
              lapply(1:J_i[i], \(j, i){
                data.frame(
                  i = i, j = j,
                  x = seq(0, 1, length.out = 100),
                  lambda = exp(beta0 + beta1 * obsCovFunctions[[i]][[j]](seq(0, 1, length.out = 100)))
                )
              }, i = i)
            })
          ),
          aes(x = x, y = lambda)
        ) +
        geom_point(
          data = bind_rows(
            lapply(1:S, \(i){
              lapply(1:J_i[i], \(j, i){
                ndet <- length(detec_list[[i]][[j]])
                data.frame(
                  i = rep(i, ndet), j = rep(j, ndet),
                  detec_time = detec_list[[i]][[j]]
                )
              }, i = i)
            })
          ), aes(y = 1, x = detec_time), color = "red"
        ) +
        theme(legend.position = "top")
    )
  }


  ## Negative log-likelihood
  negLogLik <- function(params, detec_list, obsCovFunctions, siteCov, T0 = 0, T1 = 1) {
    occu0 <- params[1]
    occu1 <- params[2]
    det0 <- params[3]
    det1 <- params[4]

    S <- length(detec_list)
    J_i <- lengths(detec_list)
    psi <- plogis(occu0 + occu1 * siteCov)
    lik_site <- rep(NA, S)

    for (i in 1:S) {
      LAMBDA_ij <- sapply(seq_len(J_i[i]), function(j) {
        res <- tryCatch(
          integrate(
            function(t) exp(det0 + det1 * obsCovFunctions[[i]][[j]](t)),
            T0, T1
          )$value,
          error = function(e) .Machine$double.xmax
        )
        res
      })
      K_i <- sum(lengths(detec_list[[i]]))
      if (K_i > 0) {
        L_i <- psi[i]
        for (j in 1:J_i[i]) {
          prodexp <- prod(exp(det0 + det1 * obsCovFunctions[[i]][[j]](detec_list[[i]][[j]])))
          if (!is.finite(prodexp)) {
            prodexp <- .Machine$double.xmax
          }
          L_i <- L_i * exp(-LAMBDA_ij[j]) * prodexp
        }
        lik_site[i] <- L_i
      } else {
        lik_site[i] <- psi[i] * prod(exp(-LAMBDA_ij)) + (1 - psi[i])
      }
    }

    -sum(log(lik_site))
  }

  ## Fit model
  fit <- optim(
    par = c(0, 0, 0, 0),
    fn = negLogLik,
    detec_list = detec_list,
    obsCovFunctions = obsCovFunctions,
    siteCov = siteCov,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  data.frame(
    rep = id_rep,
    param = c(
      "alpha0", "alpha1", "beta0", "beta1",
      "occu", "ndet", "mean_ndet_per_deploy_Z1"
    ),
    simul = c(
      alpha0, alpha1, beta0, beta1,
      mean(Z),
      sum(sapply(detec_list, \(x){
        sum(lengths(x))
      })),
      mean(sapply(detec_list[Z == 1], \(x){
        mean(lengths(x))
      }))
    ),
    estim = c(fit$par, NA, NA, NA)
  )
}, future.seed = TRUE)

bind_rows(res_fit) |>
  filter(!is.na(estim)) |>
  ggplot(aes(x = 1, y = estim)) +
  facet_wrap(~param, scales = "free") +
  geom_hline(aes(yintercept = simul), colour = "forestgreen", linewidth = 2) +
  geom_violin(alpha = .4) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


bind_rows(res_fit) |>
  filter(!is.na(estim)) |>
  ggplot(aes(x = 1, y = estim)) +
  facet_wrap(~param, scales = "free") +
  geom_hline(aes(yintercept = simul), colour = "forestgreen", linewidth = 2) +
  geom_violin(alpha = .4) +
  geom_point() +
  scale_y_log10() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


bind_rows(res_fit) |>
  filter(!is.na(estim)) |>
  mutate(bias = estim - simul, error = abs(estim - simul)) |>
  group_by(param) |>
  summarise(
    bias_min = min(bias),
    bias_mean = mean(bias),
    bias_max = max(bias),
    error_min = min(error),
    error_mean = mean(error),
    error_max = max(error)
  )
