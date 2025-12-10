# Continuous-time unmarkedFrame (time-to-each-detection)

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
    # Type checks for numeric columns
    for (col in c("begintime", "endtime")) {
      if (!inherits(object@deploymentData[[col]], "numeric")) {
        errors <- c(errors, paste0(
          "deploymentData$", col, " must be numeric"
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
    if (!inherits(object@y$obstime, "numeric")) {
      errors <- c(errors, paste0(
        "y$obstime must be numeric"
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
      if (!inherits(object@obsCovs$obsCovsMeasures$covtime, "numeric")) {
        errors <- c(errors, paste0(
          "obsCovs$obsCovsMeasures$covtime must be numeric"
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

  if (!is.null(obsCovs$obsCovsMeasures)) {
    if (!is.numeric(obsCovs$obsCovsMeasures$covtime)) {
      stop("obsCovs$obsCovsMeasures$covtime should be numeric")
    }
    if (length(setdiff(obsCovs$obsCovsMeasures$deployment, deploymentData$deployment)) > 0) {
      stop("obsCovs$obsCovsMeasures$deployment should cover deployments in deploymentData$deployment")
    }
    if (length(setdiff(deploymentData$deployment, obsCovs$obsCovsMeasures$deployment)) > 0) {
      stop("deploymentData$deployment deployments should be covered in obsCovs$obsCovsMeasures$deployment")
    }
  }
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
          x = subdf$covtime,
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
            x = subdf$covtime,
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
  if (any(!(c("deployment", "site", "begintime", "endtime") %in% colnames(deploymentData)))) {
    stop("deploymentData should have mandatory columns: deployment, site, begintime, endtime")
  }
  if (any(!(c("deployment", "obstime") %in% colnames(y)))) {
    stop("y should have mandatory columns: deployment, obstime")
  }
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


# ASSOCIATED METHODS ----

# Not relevent for continuous-time models:
# numY ;  obsToY ; obsToY<-

## Access basic data ----
# getY(umfC)

setMethod("numSites", "unmarkedFrameContinuous", function(object) {
  nrow(object@siteCovs)
})
# numSites(umfC)


setGeneric("numDeployments", function(object, ...) {
  standardGeneric("numDeployments")
})
setMethod("numDeployments", "unmarkedFrameContinuous", function(object) {
  nrow(object@deploymentData)
})
# numDeployments(umfC)

# Maximum number of deployments per site
setGeneric("maxNumDeployments", function(object, ...) {
  standardGeneric("maxNumDeployments")
})
setMethod("maxNumDeployments", "unmarkedFrameContinuous", function(object) {
  max(table(object@deploymentData$site))
})
# maxNumDeployments(umfC)


setGeneric("deploymentData", function(object, ...) {
  standardGeneric("deploymentData")
})
setMethod("deploymentData", "unmarkedFrameContinuous", function(object) {
  object@deploymentData
})
# deploymentData(umfC)

# Max number of observations per site
setMethod("obsNum", "unmarkedFrameContinuous", function(object) {
  detec_counts <- merge(
    object@deploymentData[, c("site", "deployment")],
    setNames(as.data.frame(table(object@y$deployment)), c("deployment", "nb_obs")),
    by = "deployment", all.x = TRUE
  )
  max(tapply(
    detec_counts$nb_obs,
    detec_counts$site,
    sum,
    na.rm = TRUE
  ))
})
# obsNum(umfC)


# siteCovs(umfC)


setMethod("obsCovs", "unmarkedFrameContinuous", function(object) {
  object@obsCovs
})
# obsCovs(umfC)


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
# stateCovsNames(umfC)
# stateCovsNames(umfC, detailled = FALSE)

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
# obsCovsNames(umfC)
# obsCovsNames(umfC, detailled = FALSE)


## Replace data ----
setGeneric("deploymentData<-", function(object, ...) {
  standardGeneric("deploymentData<-")
})
setReplaceMethod("deploymentData", "unmarkedFrameContinuous", function(object, value) {
  object@deploymentData <- value
  validres <- validunmarkedFrameContinuous(object)
  if (isTRUE(validres)) {
    return(object)
  } else {
    stop(paste(validres, collapse = "\n\t"))
  }
})
# deploymentData(umfC) <- NULL
# deploymentData(umfC) <- t_deploymentData[1, ]
# deploymentData(umfC) <- t_deploymentData[, 1:3]
# deploymentData(umfC) <- t_deploymentData[]




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
# obsCovs(umfC) <- NULL
# obsCovs(umfC) <- data.frame("a" = rnorm(4)) # normal error
# obsCovs(umfC) <- list(t_obsCovsMeasures, t_obsCovsFunctions)


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
# siteCovs(umfC) <- NULL
# siteCovs(umfC) <- data.frame("a" = rnorm(4)) # normal error
# siteCovs(umfC) <- list(t_obsCovsMeasures, t_obsCovsFunctions) # normal error
# siteCovs(umfC) <- t_siteCovs

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
# as(umfC, "data.frame")

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
    axis(1,
      at = seq(
        from = min(deploy$begintime),
        to = max(deploy$endtime),
        by = round((max(deploy$endtime) - min(deploy$begintime)) / ntime.ticks)
      )
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
# plot(umfC)
# plot(umfC, col.detection = "blue", show.legend = F, ntime.ticks = 2)


# show(umfC)
# print(umfC)




setMethod("summary", "unmarkedFrameContinuous", function(object, ...) {
  cat("unmarkedFrameContinuous Object\n\n")
  cat(numSites(object), "sites\n")
  cat(
    numDeployments(object), "deployments with",
    maxNumDeployments(object), "maximum per site.",
    "Total:",
    format(
      sum(object@deploymentData$endtime - object@deploymentData$begintime)
    ), "\n"
  )

  object@deploymentData$deployment_duration <- round(object@deploymentData$endtime -
    object@deploymentData$begintime, 2)
  detec_counts <- merge(
    object@deploymentData[, c("site", "deployment", "deployment_duration")],
    setNames(as.data.frame(table(object@y$deployment)), c("deployment", "nb_obs")),
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

  cat("\nObservations (y) per deployment:\n")
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
    cat("\nObservation-level covariates (estimated from interpolated functions):\n")
    t0 <- as.integer(min(object@deploymentData$begintime))
    t1 <- as.integer(max(object@deploymentData$endtime))
    covsForSummary <- lapply(names(object@obsCovs$obsCovsFunctions), function(covname) {
      unlist(lapply(object@deploymentData$deployment, function(dp) {
        covfun <- object@obsCovs$obsCovsFunctions[[covname]]
        if (!is.function(covfun)) covfun <- covfun[[dp]]
        covfun(seq(t0, t1, length.out = 10000))
      }))
    })
    names(covsForSummary) <- names(object@obsCovs$obsCovsFunctions)
    print(summary(as.data.frame(covsForSummary)))
  }
})
# summary(umfC)

# TESTS ------------------------------------------------------------------------
library(testthat)

test_that("unmarkedFrameContinuous constructor works correctly", {
  # Minimal valid input
  minimal_y <- data.frame(
    deployment = c("A1", "A1"),
    obstime = c(5, 6)
  )
  minimal_deploy <- data.frame(
    site = "A",
    deployment = c("A1", "A2"),
    begintime = c(0, 10),
    endtime = c(10, 20)
  )
  expect_no_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = NULL,
      siteCovs = NULL
    )
  })

  # Adding obsCovs: measures
  obs_meas <- data.frame(
    deployment = c("A1", "A1", "A2", "A2"),
    covtime = c(0, 10, 10, 20),
    temp = c(20, 30, 28, 19)
  )
  expect_no_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(obsCovsMeasures = obs_meas)
    )
  })

  # Adding obsCovs: function
  expect_no_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(obsCovsFunctions = list(temp = function(t) t * 2))
    )
  })
  expect_no_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(obsCovsFunctions = list(temp = list(
        A1 = function(t) 1 + t,
        A2 = function(t) 1 / t^2
      )))
    )
  })

  # y missing deployment
  expect_error(
    unmarkedFrameContinuous(
      y = data.frame(obstime = 5),
      deploymentData = minimal_deploy
    ),
    "deployment"
  )

  # y missing obstime
  expect_error(
    unmarkedFrameContinuous(
      y = data.frame(deployment = "A1"),
      deploymentData = minimal_deploy
    ),
    "obstime"
  )

  # deploymentData missing required columns
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = data.frame(
        site = "A",
        deployment = "A1",
        begintime = 0
      )
    ),
    "endtime"
  )

  # begintime >= endtime
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = data.frame(
        site = "A",
        deployment = "A1",
        begintime = 10,
        endtime = 8
      )
    ),
    "endtime"
  )

  # y$deployment not present in deploymentData
  expect_error(
    unmarkedFrameContinuous(
      y = data.frame(
        deployment = "BAD",
        obstime = 5
      ),
      deploymentData = minimal_deploy
    ),
    "deployment"
  )

  # y$obstime outside interval
  expect_error(
    unmarkedFrameContinuous(
      y = data.frame(
        deployment = "A1",
        obstime = 500
      ),
      deploymentData = minimal_deploy
    ),
    "outside"
  )

  # duplicated deployment IDs
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = data.frame(
        site = c("A", "A"),
        deployment = c("A1", "A1"),
        begintime = c(0, 1),
        endtime = c(10, 5)
      )
    ),
    "duplicated deployment"
  )

  # siteCovs missing "site"
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      siteCovs = data.frame(foo = 1:2)
    ),
    "site"
  )

  # siteCovs has duplicated site IDs
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      siteCovs = data.frame(site = c("A", "A"), x = 1:2)
    ),
    "duplicated"
  )

  # forbidden name overlap between deploymentData and siteCovs
  expect_error(
    unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      siteCovs = data.frame(site = "A", deployment = "XX")
    ),
    "common column names"
  )

  # Helper: build object but break pieces manually
  goodumf <- unmarkedFrameContinuous(minimal_y, minimal_deploy)

  ## Missing y
  expect_error(goodumf@y <- NULL)

  ## Empty y
  badumf <- goodumf
  badumf@y <- data.frame()
  expect_match(validunmarkedFrameContinuous(badumf), "missing")

  ## Missing deploymentData
  expect_error(goodumf@deploymentData <- NULL)

  ## non-numeric begintime / endtime
  expect_error(unmarkedFrameContinuous(minimal_y, data.frame(
    site = "A",
    deployment = c("A1", "A2"),
    begintime = c("0", "10"),
    endtime = c("10", "20")
  )))

  ## non-numeric obstime
  expect_error(unmarkedFrameContinuous(data.frame(
    deployment = c("A1"),
    obstime = c("5")
  ), minimal_deploy))

  ## non-numeric covtime
  expect_error({
    unmarkedFrameContinuous(minimal_y, minimal_deploy,
      obsCovs = list(data.frame(deployment = "A1", covtime = "x", x = 5))
    )
  })

  ## obsCovs: more than one df
  expect_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(
        data.frame(deployment = "A1", covtime = 1, a = 1),
        data.frame(deployment = "A1", covtime = 2, a = 2)
      )
    )
  })

  ## obsCovs: more than one function list
  expect_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(
        obsCovsFunctions = list(temp = function(t) t),
        obsCovsFunctions2 = list(temp = function(t) t)
      )
    )
  })

  ## obsCovs: functions with wrong number of args
  expect_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(list(temp = function(a, b) a + b))
    )
  })

  ## obsCovsMeasures deployments not in deploymentData
  expect_error({
    m <- unmarkedFrameContinuous(
      y = minimal_y,
      deploymentData = minimal_deploy,
      obsCovs = list(data.frame(deployment = "BAD", covtime = 1, x = 5))
    )
  })
})


test_that("unmarkedFrameContinuous associated methods works correctly", {
  # Build a clean valid object for testing
  y <- data.frame(
    deployment = c("D1", "D1", "D2"),
    obstime = c(5, 6, 15)
  )
  deploy <- data.frame(
    site = c("S1", "S1"),
    deployment = c("D1", "D2"),
    begintime = c(0, 10),
    endtime = c(10, 20)
  )
  sitec <- data.frame(
    site = "S1",
    elev = 100
  )
  obs_meas <- data.frame(
    deployment = c("D1", "D1", "D2", "D2"),
    covtime = c(0, 10, 10, 20),
    temp = c(1, 2, 3, 4)
  )
  umf <- unmarkedFrameContinuous(
    y = y,
    deploymentData = deploy,
    obsCovs = list(obsCovsMeasures = obs_meas),
    siteCovs = sitec
  )
  expect_equal(getY(umf), y)
  expect_equal(numSites(umf), 1)
  expect_equal(numDeployments(umf), 2)
  expect_equal(maxNumDeployments(umf), 2)
  expect_equal(deploymentData(umf), deploy)
  expect_equal(obsNum(umf), 3)
  expect_equal(siteCovs(umf), sitec)
  expect_true("temp" %in% names(obsCovs(umf)$obsCovsFunctions))
  expect_true("elev" %in% stateCovsNames(umf)$siteCovs)
  expect_true("elev" %in% stateCovsNames(umf, detailled = FALSE))
  expect_true("temp" %in% obsCovsNames(umf)[["obsCovs$obsCovsFunctions"]])
  expect_true("temp" %in% obsCovsNames(umf, detailled = FALSE))

  # Replace deploymentData
  new_deploy <- deploy
  new_deploy$begintime <- c(0, 11)
  expect_no_error({
    umf2 <- umf
    deploymentData(umf2) <- new_deploy
  })
  expect_error({
    umf2 <- umf
    deploymentData(umf2) <- data.frame(
      site = "S1", deployment = "D1",
      begintime = 10, endtime = 5 # begintime >= endtime
    )
  })

  # Replace obsCovs
  expect_no_error({
    umf2 <- umf
    obsCovs(umf2) <- NULL
  })
  expect_no_error({
    umf2 <- umf
    obsCovs(umf2) <- list(
      obsCovsFunctions = list(temp = function(t) t * 2)
    )
  })
  expect_error({
    umf2 <- umf
    obsCovs(umf2) <- list(
      obsCovsFunctions = list(temp = function(a, b) a + b) # Two args
    )
  })

  # Replace: siteCovs<-
  expect_no_error({
    umf2 <- umf
    siteCovs(umf2) <- NULL
  })
  expect_error({
    umf2 <- umf
    siteCovs(umf2) <- data.frame(x = 1:2)
  })

  # Coercion: from unmarkedFrameContinuous -> data.frame
  expect_no_error(df <- as(umf, "data.frame"))
  expect_true("site" %in% names(df))
  expect_true("deployment" %in% names(df))
  expect_true("nb_obs" %in% names(df))
  expect_true(any(grepl("temp.", names(df))))
  expect_error(
    # Data.frame -> unmarkedFrameContinuous: error
    as(data.frame(a = 1:3), "unmarkedFrameContinuous")
  )

  # Plot,show,summary...
  expect_no_error(plot(umf))
  expect_no_error(plot(umf, col.detection = "blue", show.legend = FALSE))
  expect_no_error(show(umf))
  expect_no_error(print(umf))
  expect_no_error(summary(umf))
})
