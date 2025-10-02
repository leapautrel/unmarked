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
  elev = c(0.8497861, 1.5632035, 0.4787604),
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

## t_obsCovsContinuous ----
t_obsCovsContinuous <- expand.grid(
  "deployment" = t_deploymentData$deployment,
  "covtime" = seq(
    lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
    lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
    by = 60 * 30
  )
)
t_obsCovsContinuous$daycov <- c(lubridate::cyclic_encoding(
  t_obsCovsContinuous$covtime,
  periods = "day", encoder = "cos"
)) + rnorm(n = length(t_obsCovsContinuous$covtime))
t_obsCovsContinuous$lincov <- NA
for (s in unique(t_deploymentData$site)) {
  tmpslope <- rnorm(1)
  for (d in unique(t_deploymentData$deployment[t_deploymentData$site == s])) {
    tmpint <- rnorm(1, sd = 5)
    tmpwhich <- which(t_obsCovsContinuous$deployment == d &
      lubridate::minute(t_obsCovsContinuous$covtime) == 0)
    t_obsCovsContinuous[tmpwhich, "lincov"] <- tmpint + c(scale(seq_along(tmpwhich))) * tmpslope +
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

## obsCovs visu ----
library(ggplot2)
ggplot(t_obsCovsContinuous, aes(x = covtime, y = daycov, colour = deployment)) +
  geom_point(alpha = .5) +
  geom_line(alpha = .5) +
  geom_smooth(method = "loess", formula = "y~x", span = 0.1, alpha = .1)
ggplot(t_obsCovsContinuous[!is.na(t_obsCovsContinuous$lincov), ], aes(x = covtime, y = lincov, colour = deployment)) +
  geom_point(alpha = .5) +
  geom_line(alpha = .5) +
  geom_smooth(method = "loess", formula = "y~x", alpha = .1)
ggplot() +
  lapply(names(t_obsCovsFunctions$funcov_dep), function(deployment) {
    geom_line(aes(
      x = seq(
        lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
        lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
        by = 60 * 30
      ),
      y = t_obsCovsFunctions$funcov_dep[[deployment]](
        as.integer(seq(
          lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
          lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
          by = 60 * 30
        ))
      ),
      colour = deployment
    ))
  }) +
  labs(x = "time", y = "t_obsCovsFunctions$funcov_dep(time)")


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

    # Only one data.frame allowed: obsCovsContinuous
    obsCovs_df_idx <- which(sapply(object@obsCovs, is.data.frame))
    if (length(obsCovs_df_idx) > 1) {
      errors <- c(
        errors,
        "obsCovs has more than one data.frame element (obsCovsContinuous)"
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

  # obsCovsContinuous if present
  if (!is.null(object@obsCovs$obsCovsContinuous)) {
    required_obsCovsCont <- c("deployment", "covtime")
    missing_obsCovsCont <- setdiff(
      required_obsCovsCont,
      names(object@obsCovs$obsCovsContinuous)
    )
    if (length(missing_obsCovsCont) > 0) {
      errors <- c(errors, paste(
        "obsCovsContinuous is missing:",
        paste(missing_obsCovsCont, collapse = ", ")
      ))
    } else {
      if (!inherits(object@obsCovs$obsCovsContinuous$covtime, "POSIXct")) {
        errors <- c(errors, paste0(
          "obsCovs$obsCovsContinuous$covtime must be POSIXct"
        ))
      }
    }
    # All deployments in obsCovsContinuous must exist in deploymentData
    local({
      mismatch <- setdiff(object@obsCovs$obsCovsContinuous$deployment, deploys)
      if (length(mismatch) > 0) {
        errors <<- c(errors, paste(
          "obsCovsContinuous contains deployment(s) not in deploymentData:",
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
unmarkedFrameContinuous <- function(y, deploymentData, siteCovs = NULL, obsCovs = NULL) {
  # Order data
  deploymentData <- deploymentData[
    order(deploymentData$site, deploymentData$begintime, deploymentData$endtime),
  ]
  y <- y[order(y$deployment, y$obstime), ]

  # List keys
  sites <- unique(deploymentData$site)
  deploys <- deploymentData$deployment


  if (!is.null(obsCovs)) {
    # Infer elements by type

    # Only one data.frame allowed: obsCovsContinuous
    obsCovs_df_idx <- which(sapply(obsCovs, is.data.frame))
    if (length(obsCovs_df_idx) == 1) {
      names(obsCovs)[obsCovs_df_idx] <- "obsCovsContinuous"
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

  umf <- new(
    Class = "unmarkedFrameContinuous",
    y = y,
    deploymentData = deploymentData,
    siteCovs = siteCovs,
    obsCovs = obsCovs
  )
  return(umf)
}

unmarkedFrameContinuous(
  y = t_y,
  deploymentData = t_deploymentData,
  siteCovs = t_siteCovs,
  obsCovs = list(t_obsCovsContinuous, t_obsCovsFunctions)
)
