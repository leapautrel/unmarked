# Fit the continuous occupancy model PP
# (Poisson Process)

# Occupancy process (single site, single species)
# z_i ~ Bernoulli(psi_i)
#
# with:
#   z_i = Occupancy state of site i
#       = 1 if the site i is occupied
#       = 0 else
#   psi_i = Occupancy probability of site i

# Detection process
# t_i | z_i = 1 ~ Poisson process (lambda(t))
# t_i | z_i = 0 ~ no detection

# with: t_i = The vector of detection times in site i ; z_i = Occupancy state of
# site i ; lambda(t) = Instantaneous detection rate depending on time ;  Poisson
# process (lambda(t)) = the time between each successive detection events
# depends on exp(lambda(t))

# Classes ----------------------------------------------------------------------
setClass("unmarkedFitContinuousPP",
  contains = "unmarkedFit"
)

# Methods ----------------------------------------------------------------------
#
# ## valeurs test
# data <- umf <- umfC
# formulas <- list(
#   psi = formula("~elev+lat"),
#   lambda = formula("~camtrap_height+funcov_cste")
# )
# psiformula <- formulas$psi
# lambdaformula <- formulas$lambda
# na.rm <- TRUE

## getDesign method ----
setMethod(
  "getDesign", "unmarkedFrameContinuous",
  function(umf, formulas, na.rm = TRUE) {
    # Retrieve informations
    M <- numSites(umf)
    sites <- siteCovs(umf)$site
    D <- numDeployments(umf)
    Dmax <- maxNumDeployments(umf)
    y <- getY(umf)
    deployData <- deploymentData(umf)
    deploys <- deployData$deployment

    deploySite <- split(deployData$deployment, deployData$site)
    kMax <- max(table(umf@y$deployment)) # Max detections/deploy

    # Potential state/obs covariates
    potObsCovs <- obsCovsNames(umf, detailled = T)

    # Clean state covariates
    sc <- cleanStateCovsContinuous(umf = umf, formulas = formulas)

    # Model matrices and offset for state submodel
    X_state <- get_model_matrix(formulas$psi, sc)
    offset_state <- get_offset(formulas$psi, sc)
    Z_state <- get_Z(formulas$psi, sc)

    # Time of each detection 3D t[i,j,k]
    obsTimes <- merge(umf@y, deployData, by = "deployment")
    obsTimes <- obsTimes[
      order(obsTimes$site, obsTimes$deployment, obsTimes$obstime),
    ]
    obsTimes$k <- ave(obsTimes$obstime, obsTimes$deployment, FUN = seq_along)
    noObs <- setdiff(deployData$deployment, unique(umf@y$deployment))
    tmp <- deployData[deployData$deployment %in% noObs, , drop = FALSE]
    tmp$obstime <- as.POSIXct(NA)
    tmp$k <- 0
    ocTimes <- rbind(obsTimes, tmp)
    ocTimes <- ocTimes[
      order(ocTimes$site, ocTimes$deployment, ocTimes$k),
    ]

    # Add functions covariates at the time of each detection
    ocFuns <- obsCovs(umf)$obsCovsFunctions
    for (ocName in names(ocFuns)) {
      ocTimes[[ocName]] <- NA_real_
      ocTimes[[ocName]] <- mapply(function(dep, t) {
        f <- if (length(ocFuns[[ocName]]) == 1) ocFuns[[ocName]] else ocFuns[[ocName]][[dep]]
        f(t)
      }, ocTimes$deployment, ocTimes$obstime)
    }

    # Model matrices and offset for detection submodel at detections
    X_det <- get_model_matrix(formulas$lambda, ocTimes)
    offset_det <- get_offset(formulas$lambda, ocTimes)
    Z_det <- get_Z(formulas$lambda, ocTimes)

    # Identify missing values in state covs
    has_na_site <- row_has_na(cbind(X_state, Z_state))
    has_na_site <- matrix(rep(has_na_site, ncol(y)), M)
    drop_sites <- row_all_na(has_na_site)

    if (na.rm & any(drop_sites)) {
      warning("Site(s) ", paste(which(drop_sites), collapse = ","),
        " dropped due to missing values",
        call. = FALSE
      )
      y <- y[!drop_sites, , drop = FALSE]
      X_state <- X_state[!drop_sites, , drop = FALSE]
      offset_state <- offset_state[!drop_sites]
      Z_state <- Z_state[!drop_sites, , drop = FALSE]
    }

    list(
      y = y,
      X_state = X_state, offset_state = offset_state, Z_state = Z_state,
      X_det = X_det, offset_det = offset_det, Z_det = Z_det,
      removed.sites = which(drop_sites)
    )
  }
)

# Function ----
occuContinuousPP <- function(
    data,
    psiformula = ~1,
    lambdaformula = ~1,
    psistarts,
    lambdastarts,
    starts,
    method = "BFGS",
    se = TRUE,
    engine = c("C", "R"),
    na.rm = TRUE,
    return.negloglik = NULL,
    n_quadra = 20,
    ...) {
  R_nll_occuContinuousPP <- function(params) {
    # Compute site occupancy probabilities
    psi <- plogis(Xpsi %*% params[psiIdx])[, 1] # occupancy prob. per site
    params_lambda <- params[lambdaIdx] # detection rate coefficients

    # Compute integral of lambda(t) for each deployment
    # LAMBDA_ij: matrix of integrated lambda per site × deployment
    LAMBDA_ij <- matrix(NA, nrow = M, ncol = Dmax)

    for (i in seq_along(sites)) {
      deploySite <- deployData$deployment[deployData$site == sites[i]]

      for (j in seq_along(deploySite)) {
        dep <- deploySite[j]
        qinfo <- detCovVal[[dep]] # precomputed quadrature info

        # Linear predictor at quadrature nodes
        eta_nodes <- as.vector(qinfo$X_nodes %*% params_lambda)
        m <- max(eta_nodes)
        lambda_nodes <- exp(eta_nodes - m)

        # Integral approximation: sum(weight * lambda(node))
        LAMBDA_ij[i, j] <- sum(qinfo$w_nodes * lambda_nodes) * exp(m)
      }
    }

    # Compute lambda at detection times
    # Store only necessary info: list of vectors per deployment
    lambda_det_list <- lapply(deploys, function(dep) {
      qinfo <- detCovVal[[dep]]
      if (!is.null(qinfo$X_detTimes)) {
        eta_det <- as.vector(qinfo$X_detTimes %*% params_lambda)
        exp(eta_det)
      } else {
        numeric(0) # no detections for this deployment
      }
    })
    names(lambda_det_list) <- deploys

    # Compute log-likelihood
    llh_vec <- numeric(M) # log-likelihood per site
    for (i in seq_along(sites)) {
      site_deps <- deployData$deployment[deployData$site == sites[i]]
      LAM <- LAMBDA_ij[i, seq_along(site_deps)] # integral λ per deployment
      lambda_dets <- lapply(site_deps, function(dep) lambda_det_list[[dep]])

      if (siteHasDet[i]) {
        # Site with at least one detection
        prod_det <- sapply(lambda_dets, function(ld) if (length(ld) == 0) 1 else prod(ld))
        llh_vec[i] <- log(psi[i] * prod(exp(-LAM) * prod_det))
      } else {
        # Site without detections
        llh_vec[i] <- log(psi[i] * prod(exp(-LAM)) + (1 - psi[i]))
      }
    }

    # Negative log-likelihood
    return(-sum(llh_vec))
  }

  # Check arguments ------------------------------------------------------------
  if (!is(data, "unmarkedFrameContinuous")) {
    stop("Data is not an unmarkedFrameContinuous object. See ?unmarkedFrameContinuous if necessary.")
  }
  stopifnot(class(psiformula) == "formula")
  stopifnot(class(lambdaformula) == "formula")
  if (!missing(psistarts)) {
    stopifnot(class(psistarts) %in% c("numeric", "double", "integer"))
  }
  if (!missing(lambdastarts)) {
    stopifnot(class(lambdastarts) %in% c("numeric", "double", "integer"))
  }
  se <- as.logical(match.arg(
    arg = as.character(se),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  na.rm <- as.logical(match.arg(
    arg = as.character(na.rm),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  engine <- match.arg(engine, c("C", "R"))

  # Do not yet manage random effects!!!
  if (has_random(psiformula) | has_random(lambdaformula)) {
    stop("occuContinuousPP does not currently handle random effects.")
  }


  # Format input data ----------------------------------------------------------
  M <- numSites(data)
  sites <- siteCovs(data)$site
  D <- numDeployments(data)
  Dmax <- maxNumDeployments(data)
  y <- getY(data)
  deployData <- deploymentData(data)
  deploys <- deployData$deployment
  potDetCov <- obsCovsNames(umfC)

  # Retrieve formlist
  formlist <- list(lambda = lambdaformula, psi = psiformula)

  # Get the design matrix (calling the getDesign method for unmarkedFrame)
  # For more informations, see: getMethod("getDesign", "unmarkedFrameOccuCOP")
  designMats <- getDesign(umf = data, formulas = formlist, na.rm = TRUE)

  # Xpsi is the fixed effects design matrix, Zpsi the random effects design matrix for occupancy
  Xpsi <- designMats$X_state
  Zpsi <- designMats$Z_state

  # Same for detection rate. Note: does not contain full information for the
  # detection covariates that vary in time
  Xlambda <- designMats$X_det
  Zlambda <- designMats$Z_det

  # Vector of the sites removed because of NAs
  sitesRemoved <- designMats$removed.sites
  removed_obs <- c()

  # Parameter indices
  ParamPsi <- colnames(Xpsi)
  NbParamPsi <- ncol(Xpsi)
  psiIdx <- 1:NbParamPsi

  ParamLambda <- colnames(Xlambda)
  NbParamLambda <- ncol(Xlambda)
  nP <- NbParamPsi + NbParamLambda
  lambdaIdx <- (NbParamPsi + 1):nP
  detCovs <- setdiff(ParamLambda, "(Intercept)")

  # Checking the starting point for optim ----
  # Check if either (psistarts AND lambdastarts) OR starts is provided
  if (!missing(psistarts) & !missing(lambdastarts)) {
    # Both psistarts and lambdastarts provided
    if (!missing(starts)) {
      if (!all(c(psistarts, lambdastarts) == starts)) {
        warning(
          "You provided psistarts, lambdastarts and starts. ",
          "Please provide either (psistarts AND lambdastarts) OR starts. ",
          "Using psistarts and lambdastarts."
        )
      }
    }
    if (length(lambdastarts) != NbParamLambda) {
      stop(
        "lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
        "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula
      )
    }
    if (length(psistarts) != NbParamPsi) {
      stop(
        "psistarts (", paste(psistarts, collapse = ", "), ") ",
        "should be of length ", NbParamPsi, " with psiformula ", psiformula
      )
    }
    starts <- c(psistarts, lambdastarts)
  } else if (!missing(starts)) {
    # starts provided
    if (length(starts) != nP) {
      stop(
        "starts (", paste(starts, collapse = ", "), ") ",
        "should be of length ", nP,
        " with psiformula ", psiformula,
        " and lambdaformula ", lambdaformula
      )
    }

    psistarts <- starts[1:NbParamPsi]
    lambdastarts <- starts[(NbParamPsi + 1):(NbParamPsi + NbParamLambda)]
  } else {
    # No arguments provided, apply default values

    if (missing(lambdastarts)) {
      # If lambda starts argument was not given:
      # 0 by default
      # so lambda = exp(0) = 1 by default
      lambdastarts <- rep(0, NbParamLambda)
    } else if (length(lambdastarts) != NbParamLambda) {
      stop(
        "lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
        "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula
      )
    }

    if (missing(psistarts)) {
      # If psi starts argument was not given
      # 0 by default
      # so psi = plogis(0) = 0.5 by default
      psistarts <- rep(0, NbParamPsi)
    } else if (length(psistarts) != NbParamPsi) {
      stop(
        "psistarts (", paste(psistarts, collapse = ", "), ") ",
        "should be of length ", NbParamPsi, " with psiformula ", psiformula
      )
    }
    starts <- c(psistarts, lambdastarts)
  }

  # Time of each detection event
  obsTimes <- merge(data@y, deployData, by = "deployment")
  obsTimes <- obsTimes[
    order(obsTimes$site, obsTimes$deployment, obsTimes$obstime),
    c("site", "deployment", "obstime")
  ]
  obsTimes$k <- ave(obsTimes$obstime, obsTimes$deployment, FUN = seq_along)

  siteHasDet <- sapply(sites, \(i){
    any(obsTimes$site == i)
  })

  # Precompute values for integral estimation
  gauss_legendre <- function(n) {
    eps <- 1e-14
    m <- floor((n + 1) / 2)
    x <- numeric(n)
    w <- numeric(n)

    for (i in 1:m) {
      z <- cos(pi * (i - 0.25) / (n + 0.5))
      z1 <- 0
      repeat {
        p1 <- 1
        p2 <- 0
        for (j in 1:n) {
          p3 <- p2
          p2 <- p1
          p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
        }
        pp <- n * (z * p1 - p2) / (z^2 - 1)
        z1 <- z
        z <- z1 - p1 / pp
        if (abs(z - z1) < eps) break
      }

      x[i] <- -z
      x[n + 1 - i] <- z
      w[i] <- 2 / ((1 - z^2) * pp^2)
      w[n + 1 - i] <- w[i]
    }

    list(nodes = x, weights = w)
  }

  # Pré-calcul des valeurs des covariables pour calculer l'intégrale
  detCovVal <- vector("list", length(deploys))
  names(detCovVal) <- deploys
  for (d in seq_along(deploys)) {
    dep <- deploys[d]
    T0 <- deployData$begintime[deployData$deployment == dep]
    T1 <- deployData$endtime[deployData$deployment == dep]

    #  Quadrature nodes
    gl <- gauss_legendre(n_quadra)
    # transform from [-1,1] -> [T0,T1]
    t_nodes <- (gl$nodes + 1) / 2 * (T1 - T0) + T0
    w_nodes <- gl$weights * (T1 - T0) / 2


    #  Build design matrix X_nodes
    n_nodes <- length(t_nodes)
    p <- 1 + length(detCovs)
    X_nodes <- matrix(0, nrow = n_nodes, ncol = p)
    colnames(X_nodes) <- c("(Intercept)", detCovs)
    X_nodes[, "(Intercept)"] <- 1

    # Fill columns
    for (lc in seq_along(detCovs)) {
      detCov <- detCovs[lc]

      if (detCov %in% potDetCov[["deploymentData"]]) {
        val <- deployData[[detCov]][deployData$deployment == dep]
        X_nodes[, detCov] <- val
      } else {
        f <- data@obsCovs$obsCovsFunctions[[detCov]]
        if (length(f) == 1) {
          X_nodes[, detCov] <- f(t_nodes)
        } else {
          X_nodes[, detCov] <- f[[dep]](t_nodes)
        }
      }
    }

    #  Detection times
    idxDet <- which(obsTimes$deployment == dep)
    det_times <- if (length(idxDet) > 0) obsTimes$obstime[idxDet] else numeric(0)

    if (length(det_times) > 0) {
      X_detTimes <- matrix(0, nrow = length(det_times), ncol = p)
      colnames(X_detTimes) <- c("(Intercept)", detCovs)
      X_detTimes[, "(Intercept)"] <- 1

      for (lc in seq_along(detCovs)) {
        detCov <- detCovs[lc]
        if (detCov %in% potDetCov[["deploymentData"]]) {
          val <- deployData[[detCov]][deployData$deployment == dep]
          X_detTimes[, detCov] <- val
        } else {
          f <- data@obsCovs$obsCovsFunctions[[detCov]]
          if (length(f) == 1) {
            X_detTimes[, detCov] <- f(det_times)
          } else {
            X_detTimes[, detCov] <- f[[dep]](det_times)
          }
        }
      }
    } else {
      X_detTimes <- NULL
    }
    detCovVal[[dep]] <- list(
      t_nodes = t_nodes, # quadrature node times (Gauss–Legendre)
      w_nodes = w_nodes, # quadrature weights for ∫ λ(t) dt
      X_nodes = X_nodes, # design matrix at node times for λ(t)
      X_detTimes = X_detTimes # design matrix at detection times (for λ(t_det))
    )
  }

  # nll function depending on engine -------------------------------------------
  if (identical(engine, "C")) {
    warning("C function for occuContinuousPP nehative log-likelihood is not yet implemented: using R")
    nll <- R_nll_occuContinuousPP
  } else if (identical(engine, "R")) {
    nll <- R_nll_occuContinuousPP
  }

  ## Run optim ----
  opt <- optim(
    starts,
    nll,
    method = method,
    hessian = se,
    ...
  )

  # Get output -----------------------------------------------------------------
  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  tmb_mod <- NULL
  fmAIC <- 2 * opt$value + 2 * nP

  # Organize effect estimates
  names(ests) <- c(ParamPsi, ParamLambda)
  psi_coef <- list(ests = ests[psiIdx], cov = as.matrix(covMat[psiIdx, psiIdx]))
  lambda_coef <- list(
    ests = ests[lambdaIdx],
    cov = as.matrix(covMat[lambdaIdx, lambdaIdx])
  )

  # No random effects
  psi_rand_info <- lambda_rand_info <- list()

  # Create unmarkedEstimates ---------------------------------------------------
  psi_uE <- unmarkedEstimate(
    name = "Occupancy probability",
    short.name = "psi",
    estimates = psi_coef$ests,
    covMat = psi_coef$cov,
    fixed = 1:NbParamPsi,
    invlink = "logistic",
    invlinkGrad = "logistic.grad",
    randomVarInfo = psi_rand_info
  )
  lambda_uE <- unmarkedEstimate(
    name = "Detection rate",
    short.name = "lambda",
    estimates = lambda_coef$ests,
    covMat = lambda_coef$cov,
    fixed = 1:NbParamLambda,
    invlink = "exp",
    invlinkGrad = "exp",
    randomVarInfo = lambda_rand_info
  )
  estimateList <- unmarkedEstimateList(list(psi = psi_uE, lambda = lambda_uE))


  # Create unmarkedFit object--------------------------------------------------
  umfit <- new(
    "unmarkedFitContinuousPP",
    fitType = "occuContinuousPP",
    call = match.call(),
    formlist = formlist,
    data = data,
    estimates = estimateList,
    sitesRemoved = sitesRemoved,
    AIC = fmAIC,
    opt = opt,
    negLogLike = opt$value,
    nllFun = nll,
    TMB = tmb_mod
  )
  return(umfit)
}


# Helper functions ----

cleanStateCovsContinuous <- function(umf, formulas) {
  errors <- c()
  depStateCovs <- intersect(
    stateCovsNames(umf, detailled = T)$deploymentData,
    all.vars(formulas$psi)
  )

  sc <- siteCovs(umf)
  deployData <- deploymentData(umf)
  sites <- sc$site
  if (length(depStateCovs) == 0) {
    return(sc)
  }

  for (cov in depStateCovs) {
    # Check if its values are unique within each site.
    # If unique: use that unique value (no message).
    cov_site_distinct <- aggregate(deployData[[cov]],
      by = list(site = deployData$site),
      FUN = function(x) length(unique(x))
    )
    non_unique_sites <- cov_site_distinct$site[cov_site_distinct$x > 1]
    if (length(non_unique_sites) == 0) {
      cov_site_val <- aggregate(deployData[[cov]],
        by = list(site = deployData$site),
        FUN = function(x) unique(x)
      )
    } else if (is.numeric(deployData[[cov]])) {
      # If non-unique and numeric: warning + keep the mean value
      warning(sprintf(
        "State covariate '%s' is not unique across deployments within site(s): %s. Using mean per site.",
        cov, paste(non_unique_sites, collapse = ", ")
      ))
      cov_site_val <- aggregate(deployData[[cov]],
        by = list(site = deployData$site),
        FUN = mean, na.rm = TRUE
      )
    } else {
      errors <- c(
        errors,
        sprintf(
          "State covariate '%s' is non-numeric and not unique within site(s): %s. Please provide siteCovs explicitly.",
          cov, paste(non_unique_sites, collapse = ", ")
        )
      )
    }
    cov_site_val <- cov_site_val[
      match(sites, cov_site_val$site), ,
      drop = FALSE
    ]
    sc[[cov]] <- cov_site_val$x
  }
  if (length(errors) > 0) {
    stop(errors)
  }
  return(sc)
}

# TESTS ----
library(testthat)

test_that("occuContinuousPP can be fitted", {
  set.seed(123)
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
  s_psi_int <- 0.5
  s_psi_covsite <- 0.4
  z <- rbinom(
    n = length(sites),
    size = 1,
    prob = plogis(s_psi_int + s_psi_covsite * df_site$covsite)
  )

  ## Detections
  s_lam_int <- -3
  s_lam_covdep <- -0.5
  s_lam_covtime <- 0.6
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
  ## umf
  umfC <- unmarkedFrameContinuous(
    y = df_y,
    deploymentData = df_dep,
    siteCovs = df_site,
    obsCovs = list(list("covtime" = fun_covtime))
  )

  expect_warning(mod1 <- occuContinuousPP(
    data = umfC,
    psiformula = ~covsite,
    lambdaformula = ~ covdep + covtime
  ), "C function for occuContinuousPP nehative log-likelihood is not yet implemented")
})

## deploymentData
ndep <- 50
df_dep <- data.frame("site" = sort(sample(letters, size = ndep, replace = T)))
df_dep$deployment <- paste0(df_dep$site, ave(df_dep$site, df_dep$site, FUN = seq_along))
df_dep$begintime <- runif(n = ndep, min = 0, max = 50)
df_dep$endtime <- df_dep$begintime + runif(n = ndep, min = 1, max = 50)
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
s_psi_int <- 0.3
s_psi_covsite <- 0.4
z <- rbinom(
  n = length(sites),
  size = 1,
  prob = plogis(s_psi_int + s_psi_covsite * df_site$covsite)
)

## Detections
s_lam_int <- -3
s_lam_covdep <- -0.5
s_lam_covtime <- 0.6
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

umf <- unmarkedFrameContinuous(
  y = df_y,
  deploymentData = df_dep,
  siteCovs = df_site,
  obsCovs = list(list("covtime" = fun_covtime))
)

# test_that("cleanStateCovsContinuous works correctly", {
#   expect_equivalent(
#     cleanStateCovsContinuous(umf,
#       formulas = list("state" = ~ covsite, "det" = ~1)
#     ),
#     siteCovs(umf)
#   )
#
#   # Works with unique values
#   expect_no_error(cleanStateCovsContinuous(umf,
#     formulas = list("state" = ~ camtrap_habitat + lon, "det" = ~1)
#   ))
#
#   # "Mean" warning with non-unique numerical values
#   expect_warning(
#     cleanStateCovsContinuous(umf,
#       formulas = list("state" = ~camtrap_height, "det" = ~1)
#     ),
#     "Using mean"
#   )
#
#   # Error with non-unique non-numerical values
#   expect_error(
#     cleanStateCovsContinuous(umf,
#       formulas = list("state" = ~camtrap_model, "det" = ~1)
#     ),
#     "non-numeric and not unique"
#   )
#
#   # Keeps site order consistent
#   deploymentData(umf) <- t_deploymentData[rev(seq_len(nrow(t_deploymentData))), ]
#   expect_equivalent(
#     cleanStateCovsContinuous(
#       umf,
#       formulas = list("state" = ~camtrap_habitat, "det" = ~1)
#     )$camtrap_habitat,
#     siteCovs(umf)$habitat
#   )
#
#   # Multiple cases
#   expect_warning(
#     expect_error(
#       cleanStateCovsContinuous(umf,
#         formulas = list(
#           "state" = ~ lon + camtrap_height + camtrap_model,
#           "det" = ~1
#         )
#       ),
#       "non-numeric and not unique",
#       fixed = TRUE
#     ),
#     "Using mean",
#     fixed = TRUE
#   )
# })
