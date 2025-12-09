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
#
# with:
#   t_i = The vector of detection times in site i
#   z_i = Occupancy state of site i
#   lambda(t) = Instantaneous detection rate depending on time
#   Poisson process (lambda(t)) =
#   the time between each successive detection events depends on exp(lambda(t))

# Methods ----------------------------------------------------------------------

## valeurs test
data <- umf <- umfC
formulas <- list(
  psi = formula("~elev+lat"),
  lambda = formula("~camtrap_height+funcov_cste")
)
psiformula <- formulas$psi
lambdaformula <- formulas$lambda
na.rm <- TRUE

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
    ...) {
  R_nll_occuContinuousPP <- function(params) {
    psi <- plogis(Xpsi %*% params[psiIdx])[, 1]
    params_lambda <- params[lambdaIdx]

    LAMBDA_ij <- matrix(NA, nrow = M, ncol = Dmax) # Integral of lambda per deployment j in site i
    lambda_ijk <- obsTimes
    lambda_ijk$lambda <- NA # lambda(t) at the time of each detection
    for (i in seq_along(sites)) {
      deploySite <- deployData$deployment[deployData$site == sites[i]]
      for (j in seq_along(deploySite)) {
        # function that returns instantaneous detection rate at x
        lambdaFun <- Vectorize(function(x) {
          vals <- numeric(length(detCovs))
          vals[1] <- 1 # Intercept
          for  (l in seq_along(detCovs)) {
            detCov <- detCovs[l]
            if (detCov %in% potDetCov[["deploymentData"]]) {
              f <- function(x) deployData[[detCov]][deployData$deployment == deploySite[j]]
            } else {
              if (length(umf@obsCovs$obsCovsFunctions[[detCov]]) == 1) {
                f <- umf@obsCovs$obsCovsFunctions[[detCov]]
              } else {
                f <- umf@obsCovs$obsCovsFunctions[[detCov]][[deploySite[j]]]
              }
            }
            vals[l + 1] <- f(x)
          }
          return(exp(sum(vals * params_lambda)))
        })

        # Detection rate integral
        T0 <- deployData$begintime[deployData$deployment == deploySite[j]]
        T1 <- deployData$endtime[deployData$deployment == deploySite[j]]
        LAMBDA_ij[i, j] <- integrate(lambdaFun, lower = T0, upper = T1)$value

        # Detection rate at detections
        idxDet <- lambda_ijk$site == sites[i] & lambda_ijk$deployment == deploySite[j]
        lambda_ijk$lambda[idxDet] <- unlist(lambdaFun(lambda_ijk$obstime[idxDet]))
      }
    }

    llh <- 0
    for (i in seq_along(sites)) {
      if (siteHasDet[i]) {
        # At least 1 detection in site i: the site is occupied
        deploySite <- deployData$deployment[deployData$site == sites[i]]
        llh <- llh + log(psi[i] * prod(
          exp(-LAMBDA_ij[i, ][!is.na(LAMBDA_ij[i, ])]) *
            sapply(seq_along(deploySite), \(j){
              idxDet <- lambda_ijk$site == sites[i] &
                lambda_ijk$deployment == deploySite[j]
              if (sum(idxDet) == 0) {
                return(1)
              } else {
                return(prod(lambda_ijk$lambda[idxDet]))
              }
            })
        ))
      } else {
        # No detection in site i: either the site is occupied and we observed no
        # detection event, or the site is not occupied
        llh <- llh + log(psi[i] * prod(
          exp(-LAMBDA_ij[i, ][!is.na(LAMBDA_ij[i, ])])
        ) + (1 - psi[i]))
      }
    }
    return(-llh)
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

  # Parameter indices
  ParamPsi <- colnames(Xpsi)
  NbParamPsi <- ncol(Xpsi)
  psiIdx <- 1:NbParamPsi

  ParamLambda <- colnames(Xlambda)
  NbParamLambda <- ncol(Xlambda)
  nP <- NbParamPsi + NbParamLambda
  lambdaIdx <- (NbParamPsi + 1):nP

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
  lambda_coef <- list(ests = ests[lambdaIdx],
                      cov = as.matrix(covMat[lambdaIdx, lambdaIdx]))
  
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
    "unmarkedFitOccuCOP",
    fitType = "occuCOP",
    call = match.call(),
    formlist = formlist,
    data = data,
    estimates = estimateList,
    sitesRemoved = sitesRemoved,
    removed_obs = removed_obs,
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

## TESTS ----
library(testthat)

t_deploymentData <- data.frame(
  deployment = c("A1", "A2", "B1", "B2", "C1"),
  site = c("A", "A", "B", "B", "C"),
  begintime = as.POSIXct("2025-09-19 00:00:00") + runif(n = 5, 0, 24 * 60 * 60),
  endtime = as.POSIXct("2025-09-26 00:00:00") + runif(n = 5, 0, 24 * 60 * 60),
  lat = c(43.442, 43.442, 43.453, 43.453, 43.426),
  lon = c(2.076, 2.076, 2.092, 2.092, 2.034),
  camtrap_model = c("Brand1", "Brand2", "Brand1", "Brand2", "Brand1"),
  camtrap_height = c(2.33, 1.44, 1.00, 1.00, 1.28),
  camtrap_habitat = c("Forest", "Forest", "Grassland", "Grassland", "City"),
  stringsAsFactors = FALSE
)


t_siteCovs <- data.frame(
  site = c("A", "B", "C"),
  habitat = c("Forest", "Grassland", "City"),
  elev = rnorm(3),
  stringsAsFactors = FALSE
)

n <- rpois(n = 1, lambda = 10)
t_y <- data.frame(
  deployment = sample(
    x = t_deploymentData$deployment,
    size = n,
    replace = TRUE,
    prob = c(0, 0, seq(1, length.out = length(t_deploymentData$deployment) - 2))
  ),
  obstime = as.POSIXct(runif(
    n = n,
    max(t_deploymentData$begintime),
    min(t_deploymentData$endtime)
  ))
)

t_obsCovsMeasures <- expand.grid(
  "deployment" = t_deploymentData$deployment,
  "covtime" = seq(
    lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
    lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
    by = 60 * 30
  )
)
t_obsCovsMeasures$daycov <- c(lubridate::cyclic_encoding(
  t_obsCovsMeasures$covtime,
  periods = "day", encoder = "cos"
)) + rnorm(n = length(t_obsCovsMeasures$covtime))
t_obsCovsMeasures$lincov <- NA
for (s in unique(t_deploymentData$site)) {
  tmpslope <- rnorm(1)
  for (d in unique(t_deploymentData$deployment[t_deploymentData$site == s])) {
    tmpint <- rnorm(1, sd = 5)
    tmpwhich <- which(t_obsCovsMeasures$deployment == d &
      lubridate::minute(t_obsCovsMeasures$covtime) == 0)
    t_obsCovsMeasures[tmpwhich, "lincov"] <- tmpint + c(scale(seq_along(tmpwhich))) * tmpslope +
      rnorm(n = length(tmpwhich))
  }
}

t_funlincov <- setNames(
  vector(mode = "list", length = nrow(t_deploymentData)),
  t_deploymentData$deployment
)
for (dp in t_deploymentData$deployment) {
  t_funlincov[[dp]] <- approxfun(
    x = as.integer(c(
      t_deploymentData$begintime[t_deploymentData$deployment == dp],
      t_deploymentData$endtime[t_deploymentData$deployment == dp]
    )),
    y = rnorm(n = 2),
    method = "linear",
    rule = 2
  )
}
t_obsCovsFunctions <- list(
  "funcov_dep" = t_funlincov,
  "funcov_cste" = approxfun(
    x = as.integer(c(
      lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
      lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour")
    )),
    y = rnorm(n = 2),
    method = "linear",
    rule = 2
  )
)


umf <- unmarkedFrameContinuous(
  y = t_y,
  deploymentData = t_deploymentData,
  siteCovs = t_siteCovs,
  obsCovs = list(t_obsCovsMeasures, t_obsCovsFunctions)
)

test_that("cleanStateCovsContinuous works correctly", {
  expect_equivalent(
    cleanStateCovsContinuous(umf,
      formulas = list("state" = ~ habitat + elev, "det" = ~1)
    ),
    siteCovs(umf)
  )

  # Works with unique values
  expect_no_error(cleanStateCovsContinuous(umf,
    formulas = list("state" = ~ camtrap_habitat + lon, "det" = ~1)
  ))

  # "Mean" warning with non-unique numerical values
  expect_warning(
    cleanStateCovsContinuous(umf,
      formulas = list("state" = ~camtrap_height, "det" = ~1)
    ),
    "Using mean"
  )

  # Error with non-unique non-numerical values
  expect_error(
    cleanStateCovsContinuous(umf,
      formulas = list("state" = ~camtrap_model, "det" = ~1)
    ),
    "non-numeric and not unique"
  )

  # Keeps site order consistent
  deploymentData(umf) <- t_deploymentData[rev(seq_len(nrow(t_deploymentData))), ]
  expect_equivalent(
    cleanStateCovsContinuous(
      umf,
      formulas = list("state" = ~camtrap_habitat, "det" = ~1)
    )$camtrap_habitat,
    siteCovs(umf)$habitat
  )

  # Multiple cases
  expect_warning(
    expect_error(
      cleanStateCovsContinuous(umf,
        formulas = list(
          "state" = ~ lon + camtrap_height + camtrap_model,
          "det" = ~1
        )
      ),
      "non-numeric and not unique",
      fixed = TRUE
    ),
    "Using mean",
    fixed = TRUE
  )
})
