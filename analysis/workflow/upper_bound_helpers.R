# Pair classes by gene and bootstrap replicate before summing their estimates.
upper_bound_pairs <- list(
  ASD = c("Combined Probands PTV", "Combined Probands Mis2"),
  DD = c("Developmental disorders PTV", "Developmental disorders Mis2")
)

validate_upper_bound_pairs <- function(models) {
  stopifnot(setequal(names(models), unlist(upper_bound_pairs)))
  for (pair in upper_bound_pairs) {
    a <- models[[pair[1]]]; b <- models[[pair[2]]]
    stopifnot(identical(rownames(a$input), rownames(b$input)),
              identical(a$draws, b$draws), identical(a$prevalence, b$prevalence),
              identical(b$scale, 1))
  }
  invisible(TRUE)
}

pool_upper_bound_draws <- function(d) {
  stopifnot(setequal(unique(d$model), unlist(upper_bound_pairs)),
            !anyDuplicated(d[c("model", "RR_cap", "replicate", "metric")]))
  result <- lapply(names(upper_bound_pairs), function(dataset) {
    pair <- upper_bound_pairs[[dataset]]
    a <- subset(d, model == pair[1])
    b <- subset(d, model == pair[2])
    keys <- c("RR_cap", "replicate", "metric")
    a <- a[do.call(order, a[keys]), ]; b <- b[do.call(order, b[keys]), ]
    rownames(a) <- rownames(b) <- NULL
    stopifnot(identical(a[keys], b[keys]))
    # Effective penetrance remains PTV-only for both disorders.
    combined <- a[ c(keys, "estimate") ]
    add <- combined$metric != "effective_penetrance"
    combined$estimate[add] <- a$estimate[add] + b$estimate[add]
    combined$dataset <- dataset
    combined
  })
  do.call(rbind, result)
}
