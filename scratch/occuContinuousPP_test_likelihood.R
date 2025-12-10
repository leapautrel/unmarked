if (interactive()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory(function() {})[1])
}

source("unmarkedFrameContinuous.R")
source("../R/All_generics.R")
source("../R/All_classes.R")
source("../R/utils.R")
source("../R/getDesign.R")
source("../R/mixedModelTools.R")
source("../R/unmarkedEstimate.R")
source("occuContinuousPP.R")

RES <- lapply(1:50, \(seed){
  ## deploymentData
  ndep <- 100
  df_dep <- data.frame("site" = sort(sample(
    sample(c(outer(letters, letters, paste0)), 100),
    size = ndep, replace = T
  )))
  df_dep$deployment <- paste0(df_dep$site, ave(df_dep$site, df_dep$site, FUN = seq_along))
  df_dep$begintime <- runif(n = ndep, min = 0, max = 100)
  df_dep$endtime <- df_dep$begintime + runif(n = ndep, min = 1, max = 100)
  df_dep$covdep <- rnorm(n = ndep)

  ## siteCovs
  sites <- unique(sort(df_dep$site))
  df_site <- data.frame(site = sites)
  df_site$covsite <- rnorm(n = length(sites))

  ## obsCovs
  tmp <- seq(floor(min(df_dep$begintime)), ceiling(max(df_dep$endtime)))
  fun_covtime <- approxfun(
    x = tmp,
    y = sin(2 * pi * tmp / max(tmp)) + rnorm(length(tmp), sd = 0.1),
    method = "linear", rule = 2
  )

  ## Occupied
  s_psi_int <- rnorm(n = 1, mean = 0.3, sd = 1)
  s_psi_covsite <- rnorm(n = 1, mean = 0, sd = 1)
  z <- rbinom(
    n = length(sites),
    size = 1,
    prob = plogis(s_psi_int + s_psi_covsite * df_site$covsite)
  )
  if (sum(z) == 0) {
    print(paste("seed:", seed, "- sum(z)=0"))
    return(NULL)
  }

  ## Detections
  s_lam_int <- rnorm(n = 1, mean = -3, sd = 1)
  s_lam_covdep <- rnorm(n = 1, mean = -0.5, sd = 1)
  s_lam_covtime <- rnorm(n = 1, mean = 0.6, sd = 1)
  df_y <- dplyr::bind_rows(
    lapply(seq_along(sites), function(i) {
      if (!z[i]) {
        return(NULL)
      }
      deploySite <- df_dep$deployment[df_dep$site == sites[i]]
      dplyr::bind_rows(
        lapply(seq_along(deploySite), function(j) {
          depIdx <- which(df_dep$deployment == deploySite[j])
          ovspg_lambda <- ceiling(max(exp(
            s_lam_int +
              s_lam_covdep * df_dep$covdep[depIdx] +
              s_lam_covtime * fun_covtime(seq(
                floor(min(df_dep$begintime)),
                ceiling(max(df_dep$endtime)),
                by = 0.1
              ))
          )))
          ovspg_n <- rpois(
            1,
            ovspg_lambda * (df_dep$endtime[depIdx] - df_dep$begintime[depIdx])
          )
          ovspg_t <- sort(runif(
            ovspg_n,
            df_dep$begintime[depIdx],
            df_dep$endtime[depIdx]
          ))
          ovspg_t_lam <- exp(
            s_lam_int +
              s_lam_covdep * df_dep$covdep[depIdx] +
              s_lam_covtime * fun_covtime(ovspg_t)
          )
          kept <- as.logical(rbinom(ovspg_n, 1, ovspg_t_lam / ovspg_lambda))
          if (sum(kept) > 0) {
            return(data.frame(
              deployment = deploySite[j],
              obstime = ovspg_t[kept]
            ))
          }
        })
      )
    })
  )
  if (nrow(df_y) < 5) {
    print(paste("seed:", seed, "- nrow(df_y)=", nrow(df_y)))
    return(NULL)
  }

  ## umf
  umfC <- unmarkedFrameContinuous(
    y = df_y,
    deploymentData = df_dep,
    siteCovs = df_site,
    obsCovs = list(list("covtime" = fun_covtime))
  )
  # summary(umfC)
  # plot(umfC)

  mod1 <- tryCatch(
    {
      occuContinuousPP(
        data = umfC,
        psiformula = ~covsite,
        lambdaformula = ~ covdep + covtime
      )
    },
    error = function(e) {
      return(NULL) # if error occurs, mod1 becomes NULL
    }
  )

  if (is.null(mod1)) {
    print(paste("seed:", seed, "is.null(mod1)"))
    return(NULL)
  }

  print(paste("seed:", seed, "ok"))

  rbind(
    data.frame(
      seed = seed,
      process = "psi",
      param = names(mod1@estimates@estimates$psi@estimates),
      simul = c(s_psi_int, s_psi_covsite),
      estim = unname(mod1@estimates@estimates$psi@estimates)
    ),
    data.frame(
      seed = seed,
      process = "lambda",
      param = names(mod1@estimates@estimates$lambda@estimates),
      simul = c(s_lam_int, s_lam_covdep, s_lam_covtime),
      estim = unname(mod1@estimates@estimates$lambda@estimates)
    )
  )
})

library(ggplot2)
library(dplyr)
library(ggh4x)
dplyr::bind_rows(RES) |>
  ggplot() +
  ggh4x::facet_nested_wrap(~ process + param, scales = "free") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(aes(x = simul, y = estim))
