# Continuous-time unmarkedFrame (time-to-each-detection)

# TEST VALUES ------------------------------------------------------------------
set.seed(123)
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
t_siteCovs <- data.frame(
  site = c("A", "B", "C"),
  habitat = c("Forest", "Grassland", "City"),
  elev = c(0.8497861, 1.5632035, 0.4787604),
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

t_obsCovsFunctions <- setNames(
  vector(mode = "list", length = nrow(t_deploymentData)),
  t_deploymentData$deployment
)
for (dp in t_deploymentData$deployment) {
  t_obsCovsFunctions[[dp]]$funlincov <- approxfun(
    x = as.integer(c(
      t_deploymentData$begintime[t_deploymentData$deployment == dp],
      t_deploymentData$endtime[t_deploymentData$deployment == dp]
    )),
    y = rnorm(n = 2),
    method = "linear",
    rule = 2
  )
}



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
  lapply(names(t_obsCovsFunctions), function(deployment) {
    geom_line(aes(
      x = seq(
        lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
        lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
        by = 60 * 30
      ),
      y = t_obsCovsFunctions[[deployment]]$funlincov(
        as.integer(seq(
          lubridate::floor_date(min(t_deploymentData$begintime), unit = "hour"),
          lubridate::ceiling_date(max(t_deploymentData$endtime), unit = "hour"),
          by = 60 * 30
        ))
      ),
      colour = deployment
    ))
  }) +
  labs(x = "time", y = "t_obsCovsFunctions$funlincov(time)")


# TEST CLASSE ------------------------------------------------------------------
validunmarkedFrameContinuous <- function(object) {
  errors <- character(0)

  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
}


setClass("unmarkedFrameGeneric",
  slots = c(y = "ANY", obsCovs = "ANY", siteCovs = "ANY")
)

setClass(
  "unmarkedFrameContinuous",
  contains = "unmarkedFrameGeneric",
  validity = validunmarkedFrameContinuous
)


# TEST CREATION FUN ------------------------------------------------------------

# Order data
t_deploymentData <- t_deploymentData[
  order(t_deploymentData$site, t_deploymentData$begintime, t_deploymentData$endtime),
]
t_y <- t_y[order(t_y$deployment, t_y$obstime), ]

# List keys
sites <- unique(t_deploymentData$site)
deploys <- t_deploymentData$deployment

# Count number of observations per site and deployment (y)
## FIXME add check if `deploy_num` column already exist
t_deploymentData$deploy_num <- as.integer(ave(
  t_deploymentData$site, t_deploymentData$site,
  FUN = seq_along
))
# deploy_nums <- unique(sort(t_deploymentData$deploy_num))
# y <- matrix(0,
#   nrow = length(sites), ncol = length(deploy_nums),
#   dimnames = list(sites, deploy_nums)
# )
# for (s in sites) {
#   tmp <- setdiff(deploy_nums, t_deploymentData$deploy_num[t_deploymentData$site == s])
#   if (length(tmp) > 0) {
#     y[s, tmp] <- NA
#   }
# }
# for (i in seq_len(nrow(t_y))) {
#   s <- t_deploymentData$site[t_deploymentData$deployment == t_y$deployment[i]]
#   d <- t_deploymentData$deploy_num[t_deploymentData$deployment == t_y$deployment[i]]
#   y[s, d] <- y[s, d] + 1
# }




new(
  Class = "unmarkedFrameContinuous",
  y = t_y,
  siteCovs = t_siteCovs,
  obsCovs = list(
    t_obsCovsContinuous,
    t_obsCovsFunctions
  )
)
