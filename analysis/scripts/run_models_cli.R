# Command-line contract for analysis/scripts/run_models.R.

run_models_help <- function() {
  paste0(
    "Options:\n",
    "  --mode <main|age|no_ces|no_overlap|test> Analysis mode (default: main)\n",
    "  --run-autism <true|false>            Run autism models in main mode\n",
    "  --run-ddd <true|false>                Run DDD models in main mode\n",
    "  --bootstraps <positive integer>       Number of bootstrap samples\n",
    "  --seed <integer>                      Random seed\n",
    "  --run-date <label>                    Shared analysis-run label\n",
    "  --autism-max-effect-size <number>     Autism component upper bound\n",
    "  --ddd-max-effect-size <number>        DDD component upper bound\n",
    "  --ces-gene-file <path>                Optional one-column CES gene table\n",
    "  --optimizer <EM|mixsqp>               Model optimizer (default: mixsqp)\n"
  )
}

parse_run_models_args <- function(args, current_date = Sys.Date()) {
  if ("--help" %in% args) return(list(help = TRUE))

  get_arg <- function(flag, default) {
    equals_match <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(equals_match) > 0L) {
      return(sub(paste0("^", flag, "="), "", equals_match[1L]))
    }
    flag_index <- match(flag, args)
    if (!is.na(flag_index)) {
      if (flag_index == length(args) ||
          startsWith(args[flag_index + 1L], "--")) {
        stop("Missing value after ", flag)
      }
      return(args[flag_index + 1L])
    }
    default
  }
  as_flag <- function(value, name) {
    normalized <- tolower(value)
    if (!normalized %in% c("true", "false")) {
      stop(name, " must be true or false.")
    }
    normalized == "true"
  }

  analysis_mode <- get_arg("--mode", "main")
  run_autism <- as_flag(get_arg("--run-autism", "true"), "--run-autism")
  run_ddd <- as_flag(get_arg("--run-ddd", "true"), "--run-ddd")
  n_bootstraps <- as.integer(get_arg("--bootstraps", "100"))
  random_seed <- as.integer(get_arg("--seed", "24312342"))
  run_date <- get_arg("--run-date", format(current_date, "%b%d_%y"))
  autism_max_effect_size <- as.numeric(
    get_arg("--autism-max-effect-size", "100")
  )
  ddd_max_effect_size <- as.numeric(get_arg("--ddd-max-effect-size", "250"))
  ces_gene_file <- get_arg("--ces-gene-file", "")
  if (!nzchar(ces_gene_file)) ces_gene_file <- NULL
  optimizer <- tolower(get_arg("--optimizer", "mixsqp"))
  optimizer <- if (optimizer == "mixsqp") {
    "mixsqp"
  } else if (optimizer == "em") {
    "EM"
  } else {
    stop("--optimizer must be EM or mixsqp.")
  }

  valid_modes <- c("main", "age", "no_ces", "no_overlap", "test")
  if (!analysis_mode %in% valid_modes) {
    stop("analysis_mode must be one of: ", paste(valid_modes, collapse = ", "))
  }
  if (!run_autism && (!run_ddd || analysis_mode != "main")) {
    stop("No analysis is selected to run.")
  }
  if (is.na(n_bootstraps) || n_bootstraps < 1L) {
    stop("--bootstraps must be a positive integer.")
  }
  if (is.na(random_seed)) stop("--seed must be an integer.")
  if (!grepl("^[A-Za-z0-9_-]+$", run_date)) {
    stop("--run-date may contain only letters, numbers, underscores, and hyphens.")
  }
  if (is.na(autism_max_effect_size) || autism_max_effect_size <= 1 ||
      is.na(ddd_max_effect_size) || ddd_max_effect_size <= 1) {
    stop("Maximum effect sizes must be numeric values greater than 1.")
  }

  list(
    help = FALSE,
    analysis_mode = analysis_mode,
    run_autism = run_autism,
    run_ddd = run_ddd,
    n_bootstraps = n_bootstraps,
    random_seed = random_seed,
    run_date = run_date,
    autism_max_effect_size = autism_max_effect_size,
    ddd_max_effect_size = ddd_max_effect_size,
    ces_gene_file = ces_gene_file,
    optimizer = optimizer
  )
}
