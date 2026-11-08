source("analysis/workflow/upper_bound_helpers.R")
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
d <- expand.grid(model = unlist(upper_bound_pairs), RR_cap = c(100, 250),
                 replicate = 0:3, metric = c("mutvar", "fraction_cases_RR5", "effective_penetrance"),
                 stringsAsFactors = FALSE)
d$estimate <- match(d$model, unlist(upper_bound_pairs)) * 10 + d$replicate
p <- pool_upper_bound_draws(d)
stopifnot(nrow(p) == 48)
for (dataset in names(upper_bound_pairs)) {
  pair <- upper_bound_pairs[[dataset]]
  for (i in which(p$dataset == dataset)) {
    x <- subset(d, model %in% pair & RR_cap == p$RR_cap[i] &
                  replicate == p$replicate[i] & metric == p$metric[i])
    expected <- if (p$metric[i] == "effective_penetrance") x$estimate[x$model == pair[1]] else sum(x$estimate)
    stopifnot(p$estimate[i] == expected)
  }
}
stopifnot(identical(p, pool_upper_bound_draws(d[nrow(d):1, ])),
          fails(pool_upper_bound_draws(d[-1, ])), fails(pool_upper_bound_draws(rbind(d, d[1, ]))))
models <- setNames(rep(list(list(input = data.frame(x = 1:3, row.names = letters[1:3]),
                               draws = matrix(1:3, 3, 100), prevalence = .01, scale = 1)), 4),
                   unlist(upper_bound_pairs))
validate_upper_bound_pairs(models)
models[[4]]$draws[1, 1] <- 2L
stopifnot(fails(validate_upper_bound_pairs(models)))
models[[4]] <- models[[3]]
models[[4]]$input <- models[[4]]$input[3:1, , drop = FALSE]
stopifnot(fails(validate_upper_bound_pairs(models)))
cat("PASS paired class sums, PTV-only penetrance, order independence, and invalid pairing guards\n")
