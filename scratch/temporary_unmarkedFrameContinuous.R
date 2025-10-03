# Continuous-time unmarkedFrame (time-to-each-detection)

# TEST VALUES ------------------------------------------------------------------
set.seed(123)

## t_deploymentData ----
t_deploymentData <- data.frame(
  deployment = c("A1", "A2", "B1", "B2", "C1"),
  site = c("A", "A", "B", "B", "C"),
  begintime = as.POSIXct("2025-09-19 00:00:00") + runif(n = 5, 0, 24 * 60 * 60),
  endtime = as.POSIXct("2025-09-26 00:00:00") + runif(n = 5, 0, 24 * 60 * 60),
  lat = c(43.442, 43.442, 43.453, 43.453, 43.426),
  lon = c(2.076, 2.076, 2.092, 2.092, 2.034),
  camtrap_model = c("Brand1", "Brand2", "Brand1", "Brand2", "Brand1"),
  camtrap_height = c(2.33, 1.44, 1.00, 1.00, 1.28),
  stringsAsFactors = FALSE
)

## t_siteCovs ----
t_siteCovs <- data.frame(
  site = c("A", "B", "C"),
  habitat = c("Forest", "Grassland", "City"),
  elev = rnorm(3),
  stringsAsFactors = FALSE
)

## y ----
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

## t_obsCovsMeasures ----
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

## t_obsCovsFunctions ----
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

# ## obsCovs visu ----
# library(ggplot2)
# ggplot(t_obsCovsMeasures, aes(x = covtime, y = daycov, colour = deployment)) +
#   geom_point(alpha = .5) +
#   geom_line(alpha = .5) +
#   geom_smooth(method = "loess", formula = "y~x", span = 0.1, alpha = .1)
# ggplot(t_obsCovsMeasures[!is.na(t_obsCovsMeasures$lincov), ], aes(x = covtime, y = lincov, colour = deployment)) +
#   geom_point(alpha = .5) +
#   geom_line(alpha = .5) +
#   geom_smooth(method = "loess", formula = "y~x", alpha = .1)
# ggplot() +
#   lapply(names(t_obsCovsFunctions$funcov_dep), function(deployment) {
#     geom_line(aes(
#       x = seq(
#         lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
#         lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
#         by = 60 * 30
#       ),
#       y = t_obsCovsFunctions$funcov_dep[[deployment]](
#         as.integer(seq(
#           lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
#           lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
#           by = 60 * 30
#         ))
#       ),
#       colour = deployment
#     ))
#   }) +
#   labs(x = "time", y = "t_obsCovsFunctions$funcov_dep(time)")


# TEST CLASSE ------------------------------------------------------------------

setClass("unmarkedFrame",
  slots = c(y = "ANY", obsCovs = "ANY", siteCovs = "ANY")
)

validunmarkedFrameContinuous <- function(object) {
  errors <- character(0)

  if (nrow(object@y) == 0) {
    return("Observation data (y) is missing.")
  }
  if (nrow(object@deploymentData) == 0) {
    return("Deployment data (deploymentData) is missing.")
  }

  # deploymentData -------------------------------------------------------------
  # must have deployment, site, begintime, endtime
  required_deploy <- c("deployment", "site", "begintime", "endtime")
  missing_deploy <- setdiff(required_deploy, names(object@deploymentData))
  if (length(missing_deploy) > 0) {
    errors <- c(errors, paste(
      "deploymentData is missing column(s):",
      paste(missing_deploy, collapse = ", ")
    ))
  } else {
    # Type checks for datetime columns
    for (col in c("begintime", "endtime")) {
      if (!inherits(object@deploymentData[[col]], "POSIXct")) {
        errors <- c(errors, paste0(
          "deploymentData$", col, " must be POSIXct"
        ))
      }
    }
  }

  deploys <- object@deploymentData$deployment
  if (any(duplicated(deploys))) {
    errors <- c(errors, paste(
      "deploymentData has duplicated deployment ids:",
      paste(deploys[duplicated(deploys)], collapse = ", ")
    ))
  }

  # y --------------------------------------------------------------------------
  # observation data = y slot (mandatory)
  required_obs <- c("deployment", "obstime")
  missing_obs <- setdiff(required_obs, names(object@y))
  if (length(missing_obs) > 0) {
    errors <- c(errors, paste(
      "y is missing column(s):",
      paste(missing_obs, collapse = ", ")
    ))
  } else {
    if (!inherits(object@y$obstime, "POSIXct")) {
      errors <- c(errors, paste0(
        "y$obstime must be POSIXct"
      ))
    }
  }

  # All deployments in y must exist in deploymentData
  local({
    mismatch <- setdiff(object@y$deployment, deploys)
    if (length(mismatch) > 0) {
      errors <<- c(errors, paste(
        "y contains deployment(s) not in deploymentData:",
        paste(sort(unique(mismatch)), collapse = ", ")
      ))
    }
  })

  # obsCovs --------------------------------------------------------------------
  if (!is.null(object@obsCovs)) {
    # Infer elements by type

    # Only one data.frame allowed: obsCovsMeasures
    obsCovs_df_idx <- which(sapply(object@obsCovs, is.data.frame))
    if (length(obsCovs_df_idx) > 1) {
      errors <- c(
        errors,
        "obsCovs has more than one data.frame element (obsCovsMeasures)"
      )
    }

    # Only one list of functions allowed: obsCovsFunctions
    obsCovs_funs_idx <- which(sapply(
      object@obsCovs,
      function(x) is.list(x) && !is.data.frame(x)
    ))
    if (length(obsCovs_funs_idx) > 1) {
      errors <- c(errors, paste(
        "obsCovs has more than one list element",
        "(should only have obsCovsFunctions, a list of functions of time)"
      ))
    }

    # Check that all elements are either data.frame or list/function
    other_idx <- setdiff(seq_along(object@obsCovs), c(obsCovs_df_idx, obsCovs_funs_idx))
    if (length(other_idx) > 0) {
      errors <- c(errors, paste0(
        "obsCovs contains element(s) of unsupported types: ",
        paste(
          paste(other_idx, names(object@obsCovs[other_idx]),
            sapply(object@obsCovs[other_idx], class),
            sep = ", "
          ),
          collapse = "; "
        )
      ))
    }
  }

  # obsCovsMeasures if present
  if (!is.null(object@obsCovs$obsCovsMeasures)) {
    required_obsCovsCont <- c("deployment", "covtime")
    missing_obsCovsCont <- setdiff(
      required_obsCovsCont,
      names(object@obsCovs$obsCovsMeasures)
    )
    if (length(missing_obsCovsCont) > 0) {
      errors <- c(errors, paste(
        "obsCovsMeasures is missing:",
        paste(missing_obsCovsCont, collapse = ", ")
      ))
    } else {
      if (!inherits(object@obsCovs$obsCovsMeasures$covtime, "POSIXct")) {
        errors <- c(errors, paste0(
          "obsCovs$obsCovsMeasures$covtime must be POSIXct"
        ))
      }
    }
    # All deployments in obsCovsMeasures must exist in deploymentData
    local({
      mismatch <- setdiff(object@obsCovs$obsCovsMeasures$deployment, deploys)
      if (length(mismatch) > 0) {
        errors <<- c(errors, paste(
          "obsCovsMeasures contains deployment(s) not in deploymentData:",
          paste(mismatch, collapse = ", ")
        ))
      }
    })
  }

  # obsCovsFunctions if present
  if (!is.null(object@obsCovs$obsCovsFunctions)) {
    namesCovsFuns <- names(object@obsCovs$obsCovsFunctions)

    for (nm in namesCovsFuns) {
      nmCovsFuns <- object@obsCovs$obsCovsFunctions[[nm]]

      # List of functions (per deployment)
      if (length(nmCovsFuns) > 1) {
        # Check coherence with deployment ids
        cov_deploys <- names(object@obsCovs$obsCovsFunctions[[nm]])
        if (!identical(sort(cov_deploys), sort(unique(deploys)))) {
          errors <- c(errors, paste0(
            "obsCovs$obsCovsFunctions$", nm,
            ": deployment ids mismatch with deploymentData$deployment"
          ))
        }

        # Check all elements are functions of time
        if (!all(sapply(nmCovsFuns, is.function))) {
          errors <- c(errors, paste0(
            "obsCovsFunctions$", nm, ": not all elements are functions"
          ))
        } else {
          if (any(sapply(nmCovsFuns, \(f) length(formals(f)) != 1))) {
            errors <- c(
              errors,
              paste0(
                "obsCovsFunctions$", nm,
                ": each function must have 1 numeric argument (time)"
              )
            )
          }
        }
      } else {
        # Single function case (constant across deployments)
        if (!is.function(nmCovsFuns)) {
          errors <- c(errors, paste0(
            "obsCovsFunctions$", nm, " is not a function"
          ))
        } else {
          if (length(formals(nmCovsFuns)) != 1) {
            errors <- c(errors, paste0(
              "obsCovsFunctions$", nm,
              ": function must have 1 numeric argument (time)"
            ))
          }
        }
      }
    }
  }

  # siteCovs -------------------------------------------------------------------
  if (!is.null(object@siteCovs)) {
    if (!("site" %in% names(object@siteCovs))) {
      errors <- c(errors, "siteCovs must contain a 'site' column")
    }
    sites <- object@siteCovs$site

    if (any(duplicated(sites))) {
      errors <- c(errors, paste(
        "siteCovs has duplicated site ids:",
        paste(sites[duplicated(sites)], collapse = ", ")
      ))
    }

    # All sites in deploymentData must exist in siteCovs if provided
    local({
      mismatch <- setdiff(object@deploymentData$site, sites)
      if (length(mismatch) > 0) {
        errors <<- c(errors, paste(
          "deploymentData has site(s) not found in siteCovs:",
          paste(mismatch, collapse = ", ")
        ))
      }
    })
  }

  # Forbidden duplicated column names ------------------------------------------
  ## Potential state covariates (in siteCovs and deploymentData)
  potStateCov <- c(colnames(object@siteCovs), colnames(object@deploymentData))
  potStateCov <- potStateCov[!potStateCov == "site"]
  if (any(duplicated(potStateCov))) {
    errors <- c(errors, paste(
      "deploymentData and siteCovs should not have common column names:",
      paste(potStateCov[duplicated(potStateCov)], collapse = ", ")
    ))
  }

  ## All detection covariates need to have been added as a function
  missingCovFun <- setdiff(
    colnames(object@obsCovs$obsCovsMeasures),
    c("deployment", "covtime", names(object@obsCovs$obsCovsFunctions))
  )
  if (length(missingCovFun)) {
    errors <- c(errors, paste(
      "obsCovs$obsCovsFunctions is missing some covariates",
      "declared in obsCovs$obsCovsMeasures:",
      paste(missingCovFun, collapse = ", ")
    ))
  }

  ## Potential detection covariates
  potDetCov <- c(
    colnames(object@deploymentData),
    names(object@obsCovs$obsCovsFunctions)
  )
  potDetCov <- potDetCov[!potDetCov == "deployment"]
  if (any(duplicated(potDetCov))) {
    errors <- c(errors, paste(
      "deploymentData and obsCovs$obsCovsFunctions",
      "should not have common names for potential detection covariates:",
      paste(potDetCov[duplicated(potDetCov)], collapse = ", ")
    ))
  }

  # Time consistency checks ----------------------------------------------------
  # deploymentData$begintime < deploymentData$endtime
  invalid_deploys <- object@deploymentData$deployment[
    object@deploymentData$begintime >= object@deploymentData$endtime
  ]
  if (length(invalid_deploys) > 0) {
    errors <- c(errors, paste(
      "deploymentData has begintime >= endtime for deployment(s):",
      paste(invalid_deploys, collapse = ", ")
    ))
  }

  # y$obstime must fall within corresponding deployment interval
  y_deploy_join <- merge(
    object@y[, c("deployment", "obstime")],
    object@deploymentData[, c("deployment", "begintime", "endtime")],
    by = "deployment", all.x = TRUE
  )
  y_time_mismatch <- y_deploy_join$deployment[
    y_deploy_join$obstime < y_deploy_join$begintime |
      y_deploy_join$obstime > y_deploy_join$endtime
  ]
  if (length(y_time_mismatch) > 0) {
    errors <- c(errors, paste(
      "y has observation times outside deployment period for deployment(s):",
      paste(unique(y_time_mismatch), collapse = ", ")
    ))
  }

  # Return errors
  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
}


setClass(
  "unmarkedFrameContinuous",
  slots = representation(
    y = "data.frame",
    deploymentData = "data.frame",
    obsCovs = "optionalList",
    siteCovs = "optionalDataFrame"
  ),
  contains = "unmarkedFrame",
  validity = validunmarkedFrameContinuous
)


# TEST CREATION FUN ------------------------------------------------------------
nameObsCovsC <- function(obsCovs) {
  if (!is.null(obsCovs)) {
    # Infer elements by type

    # Only one data.frame allowed: obsCovsMeasures
    obsCovs_df_idx <- which(sapply(obsCovs, is.data.frame))
    if (length(obsCovs_df_idx) == 1) {
      names(obsCovs)[obsCovs_df_idx] <- "obsCovsMeasures"
    }

    # Only one list of functions allowed: obsCovsFunctions
    obsCovs_funs_idx <- which(sapply(
      obsCovs,
      function(x) is.list(x) && !is.data.frame(x)
    ))
    if (length(obsCovs_funs_idx) == 1) {
      names(obsCovs)[obsCovs_funs_idx] <- "obsCovsFunctions"
    }
  }
  return(obsCovs)
}
completeObsCovsFunctions <- function(obsCovs, deploymentData) {
  deploys <- deploymentData$deployment

  namesObsCovsMeasures <- setdiff(
    colnames(obsCovs$obsCovsMeasures), c("deployment", "covtime")
  )
  namesObsCovsFunctions <- names(obsCovs$obsCovsFunctions)

  if (is.null(namesObsCovsFunctions) && is.null(namesObsCovsMeasures)) {
    return(obsCovs)
  } else {
    # Check duplicated covariates names
    dupObsCovs <- namesObsCovsMeasures[
      which(namesObsCovsMeasures %in% namesObsCovsFunctions)
    ]
    if (length(dupObsCovs) > 0) {
      namesObsCovsMeasures <- setdiff(namesObsCovsMeasures, dupObsCovs)
      warning(
        "The following covariate(s) are defined both in `obsCovs$obsCovsMeasures` ",
        "and in `obsCovs$obsCovsFunctions`: ",
        paste(dupObsCovs, collapse = ", "),
        ". Values from `obsCovs$obsCovsFunctions` will be used; ",
        "those in `obsCovs$obsCovsMeasures` will be ignored."
      )
    }

    # Create interpolations
    for (covname in namesObsCovsMeasures) {
      isCst <- all(sapply(
        split(obsCovs$obsCovsMeasures[[covname]], obsCovs$obsCovsMeasures$covtime),
        function(x) length(unique(x)) == 1
      ))

      if (isCst) {
        subdf <- obsCovs$obsCovsMeasures[
          obsCovs$obsCovsMeasures$deployment == deploys[1] &
            !is.na(obsCovs$obsCovsMeasures[[covname]]),
          c("deployment", "covtime", covname)
        ]
        covsDuringDeploy <- min(subdf$covtime) <=
          min(deploymentData$begintime) &&
          max(subdf$covtime) >=
            max(deploymentData$endtime)
        if (!covsDuringDeploy) {
          warning(
            "Covariate ", covname,
            " does not fully cover the deployment periods.",
            " First and last values will be extrapolated."
          )
        }
        obsCovs$obsCovsFunctions[[covname]] <- approxfun(
          x = as.integer(subdf$covtime),
          y = subdf[[covname]],
          method = "linear",
          rule = 2
        )
      } else {
        covsDuringAllDeploy <- character()
        obsCovs$obsCovsFunctions[[covname]] <- vector("list", length = length(deploys))
        names(obsCovs$obsCovsFunctions[[covname]]) <- deploys
        for (dp in deploys) {
          subdf <- obsCovs$obsCovsMeasures[
            obsCovs$obsCovsMeasures$deployment == dp &
              !is.na(obsCovs$obsCovsMeasures[[covname]]),
            c("deployment", "covtime", covname)
          ]
          covsDuringDeploy <- min(subdf$covtime) <=
            deploymentData[deploymentData$deployment == dp, "begintime"] &&
            max(subdf$covtime) >=
              deploymentData[deploymentData$deployment == dp, "endtime"]
          if (!covsDuringDeploy) {
            covsDuringAllDeploy <- c(covsDuringAllDeploy, dp)
          }
          obsCovs$obsCovsFunctions[[covname]][[dp]] <- approxfun(
            x = as.integer(subdf$covtime),
            y = subdf[[covname]],
            method = "linear",
            rule = 2
          )
        }
        if (length(covsDuringAllDeploy) > 0) {
          warning(
            "Covariate ", covname,
            " does not fully cover these deployment periods: ",
            paste(covsDuringAllDeploy, collapse = ", "),
            ". First and last values will be extrapolated."
          )
        }
      }
    }
  }
  return(obsCovs)
}


unmarkedFrameContinuous <- function(y, deploymentData, siteCovs = NULL, obsCovs = NULL) {
  # Order data
  deploymentData <- deploymentData[
    order(deploymentData$site, deploymentData$begintime, deploymentData$endtime),
  ]
  y <- y[order(y$deployment, y$obstime), ]

  # Create siteCovs if not provided
  if (is.null(siteCovs)) {
    siteCovs <- data.frame("site" = sort(unique(deploymentData$site)))
  }

  deploys <- deploymentData$deployment
  sites <- siteCovs$site

  # Name obsCovs elements
  obsCovs <- nameObsCovsC(obsCovs = obsCovs)

  # Create interpolation for all time-varying detection covariates
  # in object@obsCovs$obsCovsFunctions
  obsCovs <- completeObsCovsFunctions(
    obsCovs = obsCovs, deploymentData = deploymentData
  )

  # Build the unmarkedFrame object
  umf <- new(
    Class = "unmarkedFrameContinuous",
    y = y,
    deploymentData = deploymentData,
    siteCovs = siteCovs,
    obsCovs = obsCovs
  )
  return(umf)
}

umfC <- unmarkedFrameContinuous(
  y = t_y,
  deploymentData = t_deploymentData,
  siteCovs = t_siteCovs,
  obsCovs = list(t_obsCovsMeasures, t_obsCovsFunctions)
)

# ASSOCIATED METHODS ----

# For comparison with umfDiscrete : mallardUMF from unmarkedFrame help
data(mallard)
mallardUMF <- unmarkedFramePCount(
  mallard.y,
  siteCovs = mallard.site,
  obsCovs = mallard.obs
)


# Not relevent for continuous-time models:
# numY ;  obsToY ; obsToY<-

## Access basic data ----
getY(umfC)

setMethod("numSites", "unmarkedFrameContinuous", function(object) {
  nrow(object@siteCovs)
})
numSites(umfC)


setGeneric("numDeployments", function(object, ...) {
  standardGeneric("numDeployments")
})
setMethod("numDeployments", "unmarkedFrameContinuous", function(object) {
  nrow(object@deploymentData)
})
numDeployments(umfC)

# Max number of observations per site
setMethod("obsNum", "unmarkedFrameContinuous", function(object) {
  detec_counts <- merge(
    object@deploymentData[, c("site", "deployment")],
    setNames(as.data.frame(table(umfC@y$deployment)), c("deployment", "nb_obs")),
    by = "deployment", all.x = TRUE
  )
  max(tapply(
    detec_counts$nb_obs,
    detec_counts$site,
    sum,
    na.rm = TRUE
  ))
})
obsNum(umfC)


siteCovs(umfC)


setMethod("obsCovs", "unmarkedFrameContinuous", function(object) {
  object@obsCovs
})
obsCovs(umfC)


## Access potential covariates ----
# Extract potential state covs names from an unmarkedFrameContinuous
# Method in unmarkedFrame.R
setGeneric("stateCovsNames", function(object, ...) {
  standardGeneric("stateCovsNames")
})
setMethod("stateCovsNames", "unmarkedFrameContinuous", function(object, detailled = TRUE) {
  stateCovsNames <- vector(mode = "list", length = 2)
  names(stateCovsNames) <- c(
    "deploymentData", "siteCovs"
  )
  stateCovsNames[["deploymentData"]] <- setdiff(
    colnames(object@deploymentData),
    c("deployment", "site", "begintime", "endtime")
  )
  stateCovsNames[["siteCovs"]] <- setdiff(
    colnames(object@siteCovs),
    c("site")
  )
  if (detailled) {
    return(stateCovsNames)
  } else {
    return(unname(unlist(stateCovsNames)))
  }
})
stateCovsNames(umfC)
stateCovsNames(umfC, detailled = FALSE)

# Extract potential state covs names from an unmarkedFrameContinuous
# Method in unmarkedFrame.R
setGeneric("obsCovsNames", function(object, ...) {
  standardGeneric("obsCovsNames")
})
setMethod("obsCovsNames", "unmarkedFrameContinuous", function(object, detailled = TRUE) {
  obsCovsNames <- vector(mode = "list", length = 2)
  names(obsCovsNames) <- c(
    "deploymentData",
    "obsCovs$obsCovsFunctions"
  )
  obsCovsNames[["deploymentData"]] <- setdiff(
    colnames(object@deploymentData),
    c("deployment", "site", "begintime", "endtime")
  )

  if (!is.null(object@obsCovs$obsCovsFunctions)) {
    obsCovsNames[["obsCovs$obsCovsFunctions"]] <- names(object@obsCovs$obsCovsFunctions)
  }

  if (detailled) {
    return(obsCovsNames)
  } else {
    return(unname(unlist(obsCovsNames)))
  }
})
obsCovsNames(umfC)
obsCovsNames(umfC, detailled = FALSE)


## Replace covariates data ----
setReplaceMethod("obsCovs", "unmarkedFrameContinuous", function(object, value) {
  object@obsCovs <- completeObsCovsFunctions(
    obsCovs = nameObsCovsC(value),
    deploymentData = object@deploymentData
  )
  validres <- validunmarkedFrameContinuous(object)
  if (isTRUE(validres)) {
    return(object)
  } else {
    stop(paste(validres, collapse = "\n\t"))
  }
})
obsCovs(umfC) <- NULL
# obsCovs(umfC) <- data.frame("a" = rnorm(4)) # normal error
obsCovs(umfC) <- list(t_obsCovsMeasures, t_obsCovsFunctions)


setReplaceMethod("siteCovs", "unmarkedFrameContinuous", function(object, value) {
  if (is.null(value)) {
    value <- data.frame(
      "site" = sort(unique(object@deploymentData$site))
    )
  }
  object@siteCovs <- value
  validres <- validunmarkedFrameContinuous(object)
  if (isTRUE(validres)) {
    return(object)
  } else {
    stop(paste(validres, collapse = "\n\t"))
  }
})
siteCovs(umfC) <- NULL
# siteCovs(umfC) <- data.frame("a" = rnorm(4)) # normal error
# siteCovs(umfC) <- list(t_obsCovsMeasures, t_obsCovsFunctions) # normal error
siteCovs(umfC) <- t_siteCovs

## Coercion ----
setAs("data.frame", "unmarkedFrameContinuous", function(from) {
  stop("A data.frame can not be converted to an unmarkedFrameContinuous object.")
})
# as(data.frame("a" = rnorm(3)), "unmarkedFrameContinuous") # normal error


setAs("unmarkedFrameContinuous", "data.frame", function(from) {
  umfC <- from

  ## left join site with deployments
  umfCdf <- merge(umfC@siteCovs, umfC@deploymentData,
    by = "site", all.x = TRUE
  )

  ## Add count of detections per deployment
  detec_counts <- as.data.frame(table(umfC@y$deployment))
  colnames(detec_counts) <- c("deployment", "nb_obs")
  umfCdf <- merge(umfCdf, detec_counts, by = "deployment", all.x = TRUE)
  umfCdf$nb_obs[is.na(umfCdf$nb_obs)] <- 0

  # Add obsCovs summary (estimated with 10000 points for each cov/deployment)
  smry_suffix <- c("min", "1stQu", "median", "mean", "3rdQu", "max")
  if (!is.null(umfC@obsCovs)) {
    for (covname in names(umfC@obsCovs$obsCovsFunctions)) {
      for (dp in umfC@deploymentData$deployment) {
        whichdp <- which(umfCdf$deployment == dp)
        if (is.function(umfC@obsCovs$obsCovsFunctions[[covname]])) {
          covfun <- umfC@obsCovs$obsCovsFunctions[[covname]]
        } else {
          covfun <- umfC@obsCovs$obsCovsFunctions[[covname]][[dp]]
        }
        t0 <- as.integer(umfC@deploymentData$begintime[
          umfC@deploymentData$deployment == dp
        ])
        t1 <- as.integer(umfC@deploymentData$endtime[
          umfC@deploymentData$deployment == dp
        ])
        smry <- summary(covfun(seq(t0, t1, length.out = 10000)))
        for (i in seq_along(smry_suffix)) {
          umfCdf[whichdp, paste0(covname, ".", smry_suffix[i])] <- smry[i]
        }
      }
    }
  }

  # Reorder columns
  umfCdf <- umfCdf[, c(
    "site", "deployment", "nb_obs",
    setdiff(names(umfCdf), c("site", "deployment", "nb_obs"))
  )]
  return(umfCdf)
})
as(mallardUMF, "data.frame")
as(umfC, "data.frame")

## Display: plot/show/summary ----

setMethod(
  "plot", c(x = "unmarkedFrameContinuous"),
  function(x,
           col.deployment = adjustcolor("grey20", alpha.f = 0.15),
           col.detection = "red",
           lwd.deployment = 8,
           pch.detection = 4,
           show.legend = TRUE,
           time.format = "%Y-%m-%d",
           ntime.ticks = 10,
           ...) {
    ## Extract object info
    sites <- x@siteCovs$site
    deploy <- x@deploymentData
    y <- x@y

    site_pos <- setNames(seq_along(sites), sites)

    ## Empty plot
    plot(NA, NA,
      xlim = range(c(deploy$begintime, deploy$endtime, y$obstime)),
      ylim = c(0.5, length(sites) + 0.5),
      yaxt = "n", xaxt = "n",
      xlab = "", ylab = ""
    )

    ## Site labels on y
    axis(2, at = seq_along(sites), labels = sites, las = 1)

    ## Time axis
    axis.POSIXct(1,
      at = as.POSIXct(seq(
        from = as.Date(min(deploy$begintime)),
        to = as.Date(max(deploy$endtime)),
        by = round(as.numeric(
          difftime(max(deploy$endtime), min(deploy$begintime), units = "days")
        ) / ntime.ticks)
      )),
      format = time.format
    )


    ## Deployments
    for (i in seq_len(nrow(deploy))) {
      y0 <- site_pos[deploy$site[i]]
      segments(deploy$begintime[i], y0,
        deploy$endtime[i], y0,
        col = col.deployment, lwd = lwd.deployment
      )
    }

    ## Detections
    for (i in seq_len(nrow(y))) {
      dep <- y$deployment[i]
      site_i <- deploy$site[deploy$deployment == dep][1]
      y0 <- site_pos[site_i]
      points(y$obstime[i], y0,
        pch = pch.detection, col = col.detection,
        cex = 1.2, lwd = 1.5
      )
    }

    ## Legend
    if (show.legend) {
      legend("top",
        horiz = TRUE,
        legend = c("Deployment", "Detection"),
        col = c(col.deployment, col.detection),
        lwd = c(lwd.deployment, NA),
        pch = c(NA, pch.detection),
        pt.cex = 1.2, bty = "n"
      )
    }
  }
)
plot(umfC)
plot(umfC, col.detection = "blue", show.legend = F)


show(umfC)
print(umfC)




setMethod("summary", "unmarkedFrameContinuous", function(object, ...) {
  cat("unmarkedFrameContinuous Object\n\n")
  cat(numSites(object), "sites\n")
  cat(
    numDeployments(object), "deployments:",
    format(
      sum(object@deploymentData$endtime - object@deploymentData$begintime)
    ), "\n"
  )

  object@deploymentData$deployment_duration <- round(object@deploymentData$endtime -
    object@deploymentData$begintime, 2)
  detec_counts <- merge(
    object@deploymentData[, c("site", "deployment", "deployment_duration")],
    setNames(as.data.frame(table(umfC@y$deployment)), c("deployment", "nb_obs")),
    by = "deployment", all.x = TRUE
  )
  detec_counts$nb_obs[is.na(detec_counts$nb_obs)] <- 0

  detec_counts_site <- tapply(
    detec_counts$nb_obs,
    detec_counts$site,
    sum,
    na.rm = TRUE
  )
  cat(
    "Maximum number of observations per site:",
    max(detec_counts_site), "\n"
  )
  cat(
    "Mean number of observations per site:",
    round(mean(detec_counts_site), 2), "\n"
  )
  cat(
    "Sites with at least one detection:",
    sum(detec_counts_site > 0), "\n"
  )

  cat("y observations per deployment:")
  print(detec_counts, exclude = NULL)


  if (!is.null(object@siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(
      object@siteCovs[, -which(names(object@siteCovs) == "site"),
        drop = FALSE
      ]
    ))
  }

  if (ncol(object@deploymentData) > 4) {
    cat("\nDeployment-level covariates:\n")
    print(summary(
      object@deploymentData[,
        -which(names(object@deploymentData) %in%
          c("deployment", "site", "begintime", "endtime")),
        drop = FALSE
      ]
    ))
  }

  if (!is.null(object@obsCovs)) {
    cat("\nObservation-level covariates (estimate from interpolated functions):\n")
    t0 <- as.integer(min(umfC@deploymentData$begintime))
    t1 <- as.integer(max(umfC@deploymentData$endtime))
    covsForSummary <- lapply(names(umfC@obsCovs$obsCovsFunctions), function(covname) {
      unlist(lapply(umfC@deploymentData$deployment, function(dp) {
        covfun <- umfC@obsCovs$obsCovsFunctions[[covname]]
        if (!is.function(covfun)) covfun <- covfun[[dp]]
        covfun(seq(t0, t1, length.out = 10000))
      }))
    })
    names(covsForSummary) <- names(umfC@obsCovs$obsCovsFunctions)
    print(summary(as.data.frame(covsForSummary)))
  }
})
summary(umfC)
