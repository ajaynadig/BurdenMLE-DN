# Validated model-artifact selection for the manuscript analysis workflow.

MODEL_MANIFEST_SCHEMA_VERSION <- 1L

model_artifact_specs <- list(
  autism = list(
    model_object = "BurdenMLE_DN_models_autism",
    required_objects = c("BurdenMLE_DN_models_autism", "autism_data")
  ),
  ddd = list(
    model_object = "BurdenMLE_DN_models_DDD",
    required_objects = c("BurdenMLE_DN_models_DDD", "kaplanis_data")
  ),
  age = list(
    model_object = "BurdenMLE_DN_models_autism_AgePhenos",
    required_objects = c(
      "BurdenMLE_DN_models_autism_AgePhenos", "autism_data_AgePhenos",
      "autism_counts_AgePhenos", "age_phenotypes_spark"
    )
  ),
  no_ces = list(
    model_object = "BurdenMLE_DN_models_autism_noCES",
    required_objects = c(
      "BurdenMLE_DN_models_autism_noCES", "autism_data", "ces_gene_ids"
    )
  ),
  no_overlap = list(
    model_object = "BurdenMLE_DN_models_autism_NoKaplanis",
    required_objects = c(
      "BurdenMLE_DN_models_autism_NoKaplanis", "autism_data_NoKaplanis",
      "autism_counts_NoKaplanis"
    )
  ),
  test = list(
    model_object = "BurdenMLE_DN_models_autism_test",
    required_objects = c("BurdenMLE_DN_models_autism_test", "autism_data")
  )
)

artifact_spec <- function(artifact) {
  spec <- model_artifact_specs[[artifact]]
  if (is.null(spec)) {
    stop("Unknown model artifact '", artifact, "'.", call. = FALSE)
  }
  spec
}

is_missing_path <- function(path) {
  is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)
}

validate_character_scalar <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(name, " must be one nonempty string.", call. = FALSE)
  }
  value
}

collect_artifact_fits <- function(models, path = character()) {
  if (inherits(models, "BurdenMLEDN_fit")) {
    return(setNames(list(models), paste(path, collapse = "/")))
  }
  if (length(models) == 1L && is.atomic(models) && is.na(models)) {
    return(list())
  }
  if (!is.list(models)) {
    stop(
      "Fitted-model container has an unexpected value at ",
      paste(path, collapse = "/"), ".", call. = FALSE
    )
  }
  result <- list()
  for (i in seq_along(models)) {
    result <- c(result, collect_artifact_fits(models[[i]], c(path, i)))
  }
  result
}

validate_fit_status <- function(status, fit_id, require_convergence = TRUE) {
  required <- c("optimizer", "code", "converged", "usable")
  if (!is.list(status) || any(!required %in% names(status))) {
    stop("Fit ", fit_id, " has missing or invalid fit_status.", call. = FALSE)
  }
  if (!is.character(status$optimizer) || length(status$optimizer) != 1L ||
      is.na(status$optimizer) || !nzchar(status$optimizer) ||
      !is.character(status$code) || length(status$code) != 1L ||
      is.na(status$code) || !nzchar(status$code) ||
      !is.logical(status$converged) || length(status$converged) != 1L ||
      is.na(status$converged) ||
      !is.logical(status$usable) || length(status$usable) != 1L ||
      is.na(status$usable)) {
    stop("Fit ", fit_id, " has malformed fit_status fields.", call. = FALSE)
  }
  if (!isTRUE(status$usable)) {
    stop("Fit ", fit_id, " is unusable (", status$code, ").", call. = FALSE)
  }
  if (require_convergence && !isTRUE(status$converged)) {
    stop("Fit ", fit_id, " did not converge (", status$code, ").", call. = FALSE)
  }
}

inspect_uncertainty_stage <- function(output, fit_id, stage) {
  statuses <- output$fit_status
  if (!is.list(statuses) || !length(statuses)) {
    stop("Fit ", fit_id, " has invalid ", stage, " status.", call. = FALSE)
  }
  status_ids <- names(statuses)
  if (is.null(status_ids)) status_ids <- seq_along(statuses)
  for (i in seq_along(statuses)) {
    validate_fit_status(
      statuses[[i]], paste0(fit_id, "/", stage, "/", status_ids[i]),
      require_convergence = FALSE
    )
  }
  reliable <- all(vapply(
    statuses, function(status) isTRUE(status$converged), logical(1)
  ))
  if (!is.logical(output$reliable) || length(output$reliable) != 1L ||
      is.na(output$reliable) || !identical(output$reliable, reliable)) {
    stop("Fit ", fit_id, " has inconsistent ", stage, " reliability.",
         call. = FALSE)
  }
  list(reliable = reliable, replicates = length(statuses))
}

inspect_artifact_environment <- function(envir, artifact, require_status = TRUE) {
  spec <- artifact_spec(artifact)
  missing_objects <- spec$required_objects[
    !vapply(spec$required_objects, exists, logical(1), envir = envir,
            inherits = FALSE)
  ]
  if (length(missing_objects)) {
    stop(
      "Artifact '", artifact, "' is missing required objects: ",
      paste(missing_objects, collapse = ", "), ".", call. = FALSE
    )
  }

  models <- get(spec$model_object, envir = envir, inherits = FALSE)
  if (!is.list(models)) {
    stop("Artifact '", artifact, "' has an invalid model container.",
         call. = FALSE)
  }
  if (!require_status) {
    return(list(required_objects = spec$required_objects))
  }

  fits <- collect_artifact_fits(models)
  if (!length(fits)) {
    stop("Artifact '", artifact, "' contains no fitted models.", call. = FALSE)
  }
  fit_ids <- names(fits)
  gene_ids <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    required_fit_fields <- c(
      "component_endpoints", "delta", "features", "conditional_likelihood",
      "conditional_log_likelihood", "likelihood_log_scale", "grid_size",
      "fit_status", "uncertainty_reliable"
    )
    if (any(!required_fit_fields %in% names(fit))) {
      stop("Fit ", fit_ids[i], " is missing required modern fit fields.",
           call. = FALSE)
    }
    validate_fit_status(fit$fit_status, fit_ids[i])
    genes <- rownames(fit$features)
    if (!is.matrix(fit$features) || !is.numeric(fit$features) ||
        anyNA(fit$features) || any(!is.finite(fit$features)) ||
        any(!fit$features %in% c(0, 1)) ||
        any(rowSums(fit$features) != 1) || ncol(fit$features) < 1L ||
        is.null(colnames(fit$features)) || anyNA(colnames(fit$features)) ||
        any(!nzchar(colnames(fit$features))) ||
        anyDuplicated(colnames(fit$features)) ||
        is.null(genes) || length(genes) < 2L || anyNA(genes) ||
        any(!nzchar(genes)) || anyDuplicated(genes)) {
      stop("Fit ", fit_ids[i], " has an invalid feature matrix.",
           call. = FALSE)
    }
    endpoints <- fit$component_endpoints
    if (!is.numeric(endpoints) || !is.null(dim(endpoints)) ||
        length(endpoints) < 2L || anyNA(endpoints) ||
        any(!is.finite(endpoints)) || any(diff(endpoints) <= 0)) {
      stop("Fit ", fit_ids[i], " has invalid component endpoints.",
           call. = FALSE)
    }
    delta <- fit$delta
    expected_delta_dim <- c(ncol(fit$features), length(endpoints))
    if (!is.matrix(delta) || !is.numeric(delta) ||
        !identical(dim(delta), expected_delta_dim) || anyNA(delta) ||
        any(!is.finite(delta)) || any(delta < 0) ||
        any(abs(rowSums(delta) - 1) > 1e-8) ||
        !identical(rownames(delta), colnames(fit$features))) {
      stop("Fit ", fit_ids[i], " has invalid mixture weights.",
           call. = FALSE)
    }
    if (!identical(fit$fit_status$weights, delta)) {
      stop("Fit ", fit_ids[i],
           " has mixture weights inconsistent with fit_status.",
           call. = FALSE)
    }
    likelihood_fields <- c("conditional_likelihood", "conditional_log_likelihood")
    for (field in likelihood_fields) {
      value <- fit[[field]]
      if (!is.matrix(value) || !is.numeric(value) ||
          !identical(dim(value), c(length(genes), length(endpoints))) ||
          !identical(rownames(value), genes)) {
        stop("Fit ", fit_ids[i], " has inconsistent gene identities in ",
             field, ".", call. = FALSE)
      }
      if (identical(field, "conditional_likelihood") &&
          (anyNA(value) || any(!is.finite(value)) || any(value < 0))) {
        stop("Fit ", fit_ids[i], " has invalid conditional likelihoods.",
             call. = FALSE)
      }
      if (identical(field, "conditional_log_likelihood") &&
          (anyNA(value) || any(value == Inf))) {
        stop("Fit ", fit_ids[i], " has invalid conditional log likelihoods.",
             call. = FALSE)
      }
    }
    if (!is.numeric(fit$likelihood_log_scale) ||
        length(fit$likelihood_log_scale) != length(genes) ||
        anyNA(fit$likelihood_log_scale) ||
        any(!is.finite(fit$likelihood_log_scale)) ||
        !identical(names(fit$likelihood_log_scale), genes)) {
      stop("Fit ", fit_ids[i],
           " has inconsistent gene identities in likelihood_log_scale.",
           call. = FALSE)
    }
    genes
  })
  if (!all(vapply(gene_ids[-1L], identical, logical(1), gene_ids[[1L]]))) {
    stop("Fitted models do not share one exact gene universe and order.",
         call. = FALSE)
  }

  records <- do.call(rbind, lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    reliable <- fit$uncertainty_reliable
    if (is.null(reliable)) reliable <- NA
    if (!is.logical(reliable) || length(reliable) != 1L) {
      stop("Fit ", fit_ids[i], " has invalid uncertainty_reliable.",
           call. = FALSE)
    }
    stages <- list()
    if (!is.null(fit$bootstrap_output)) {
      stages$bootstrap <- inspect_uncertainty_stage(
        fit$bootstrap_output, fit_ids[i], "bootstrap"
      )
    }
    if (!is.null(fit$null_output)) {
      stages$null <- inspect_uncertainty_stage(
        fit$null_output, fit_ids[i], "null"
      )
    }
    expected_reliable <- if (!length(stages)) {
      NA
    } else {
      all(vapply(stages, `[[`, logical(1), "reliable"))
    }
    if (!identical(reliable, expected_reliable)) {
      stop("Fit ", fit_ids[i],
           " has inconsistent aggregate uncertainty reliability.",
           call. = FALSE)
    }
    bootstrap_count <- if (is.null(stages$bootstrap)) {
      0L
    } else {
      stages$bootstrap$replicates
    }
    if (!is.numeric(fit$grid_size) || length(fit$grid_size) != 1L ||
        is.na(fit$grid_size) || !is.finite(fit$grid_size) ||
        fit$grid_size < 1 || fit$grid_size != as.integer(fit$grid_size)) {
      stop("Fit ", fit_ids[i], " has invalid grid_size.", call. = FALSE)
    }
    maximum_effect_size <- exp(max(fit$component_endpoints))
    if (!is.finite(maximum_effect_size)) {
      stop("Fit ", fit_ids[i], " has a nonfinite maximum effect size.",
           call. = FALSE)
    }
    data.frame(
      fit_id = fit_ids[i],
      optimizer = fit$fit_status$optimizer,
      status_code = fit$fit_status$code,
      converged = fit$fit_status$converged,
      usable = fit$fit_status$usable,
      uncertainty_reliable = reliable,
      bootstrap_replicates = as.integer(bootstrap_count),
      grid_size = as.integer(fit$grid_size),
      maximum_effect_size = maximum_effect_size,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(records) <- NULL
  list(
    required_objects = spec$required_objects,
    gene_ids = gene_ids[[1L]],
    fit_records = records
  )
}

inspect_model_artifact <- function(path, artifact, require_status = TRUE) {
  validate_character_scalar(path, "Artifact path")
  if (!file.exists(path)) stop("Model artifact does not exist: ", path,
                               call. = FALSE)
  normalized <- normalizePath(path, mustWork = TRUE)
  loaded <- new.env(parent = emptyenv())
  load(normalized, envir = loaded)
  inspection <- inspect_artifact_environment(
    loaded, artifact, require_status = require_status
  )
  inspection$path <- normalized
  inspection$environment <- loaded
  inspection
}

repository_state <- function(repo_dir) {
  commit <- tryCatch(
    system2("git", c("-C", repo_dir, "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  dirty <- tryCatch(
    system2("git", c("-C", repo_dir, "status", "--porcelain",
                     "--untracked-files=no"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(commit) != 1L || !nzchar(commit)) {
    stop("Could not determine the repository commit.", call. = FALSE)
  }
  list(commit = unname(commit), tracked_dirty = length(dirty) > 0L)
}

write_model_manifest <- function(path, run_id, mode, repo_dir, settings,
                                 artifacts) {
  validate_character_scalar(path, "Manifest path")
  validate_character_scalar(run_id, "run_id")
  validate_character_scalar(mode, "mode")
  if (!is.list(settings) || is.null(names(settings)) || any(!nzchar(names(settings)))) {
    stop("settings must be a named list.", call. = FALSE)
  }
  if (!is.character(artifacts) || !length(artifacts) ||
      is.null(names(artifacts)) || any(!nzchar(names(artifacts))) ||
      anyDuplicated(names(artifacts))) {
    stop("artifacts must be a uniquely named character vector.", call. = FALSE)
  }
  entries <- lapply(seq_along(artifacts), function(i) {
    inspected <- inspect_model_artifact(artifacts[[i]], names(artifacts)[i])
    inspected$environment <- NULL
    inspected
  })
  names(entries) <- names(artifacts)
  manifest <- list(
    schema_version = MODEL_MANIFEST_SCHEMA_VERSION,
    run_id = run_id,
    mode = mode,
    repository = repository_state(repo_dir),
    settings = settings,
    artifacts = entries
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(manifest, path)
  invisible(read_model_manifest(path))
  normalizePath(path, mustWork = TRUE)
}

read_model_manifest <- function(path) {
  validate_character_scalar(path, "Manifest path")
  if (!file.exists(path)) stop("Model manifest does not exist: ", path,
                               call. = FALSE)
  manifest <- readRDS(path)
  required <- c("schema_version", "run_id", "mode", "repository", "settings",
                "artifacts")
  if (!is.list(manifest) || any(!required %in% names(manifest))) {
    stop("Model manifest is malformed.", call. = FALSE)
  }
  if (!identical(manifest$schema_version, MODEL_MANIFEST_SCHEMA_VERSION)) {
    stop("Unsupported model manifest schema version: ",
         paste(manifest$schema_version, collapse = ", "), ".", call. = FALSE)
  }
  validate_character_scalar(manifest$run_id, "Manifest run_id")
  validate_character_scalar(manifest$mode, "Manifest mode")
  allowed_artifacts <- list(
    main = c("autism", "ddd"), age = "age", no_ces = "no_ces",
    no_overlap = "no_overlap", test = "test"
  )
  if (!manifest$mode %in% names(allowed_artifacts)) {
    stop("Manifest mode is unsupported: ", manifest$mode, ".", call. = FALSE)
  }
  repository <- manifest$repository
  if (!is.list(repository) ||
      !all(c("commit", "tracked_dirty") %in% names(repository))) {
    stop("Manifest repository metadata are malformed.", call. = FALSE)
  }
  validate_character_scalar(repository$commit, "Manifest repository commit")
  if (!is.logical(repository$tracked_dirty) ||
      length(repository$tracked_dirty) != 1L || is.na(repository$tracked_dirty)) {
    stop("Manifest tracked-dirty status is malformed.", call. = FALSE)
  }
  settings <- manifest$settings
  required_settings <- c(
    "optimizer", "seed", "grid_size", "bootstrap_count", "max_effect_size"
  )
  if (!is.list(settings) || any(!required_settings %in% names(settings))) {
    stop("Manifest settings are malformed.", call. = FALSE)
  }
  validate_character_scalar(settings$optimizer, "Manifest optimizer")
  if (!settings$optimizer %in% c("mixsqp", "EM")) {
    stop("Manifest optimizer is unsupported.", call. = FALSE)
  }
  integer_setting <- function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < minimum || value != as.integer(value)) {
      stop(name, " is malformed.", call. = FALSE)
    }
  }
  integer_setting(settings$seed, "Manifest seed", -.Machine$integer.max)
  integer_setting(settings$grid_size, "Manifest grid size", 1L)
  integer_setting(settings$bootstrap_count, "Manifest bootstrap count", 0L)
  if (!is.list(settings$max_effect_size) ||
      is.null(names(settings$max_effect_size)) ||
      any(!nzchar(names(settings$max_effect_size))) ||
      any(!vapply(settings$max_effect_size, function(value) {
        is.numeric(value) && length(value) == 1L && !is.na(value) &&
          is.finite(value) && value > 1
      }, logical(1)))) {
    stop("Manifest maximum-effect settings are malformed.", call. = FALSE)
  }
  if (!is.list(manifest$artifacts) || is.null(names(manifest$artifacts)) ||
      any(!nzchar(names(manifest$artifacts))) || anyDuplicated(names(manifest$artifacts))) {
    stop("Model manifest artifacts are malformed.", call. = FALSE)
  }
  if (!length(manifest$artifacts) ||
      any(!names(manifest$artifacts) %in% allowed_artifacts[[manifest$mode]])) {
    stop("Manifest artifacts do not match its analysis mode.", call. = FALSE)
  }
  manifest
}

copy_artifact_objects <- function(from, to, object_names) {
  for (object_name in object_names) {
    assign(object_name, get(object_name, envir = from, inherits = FALSE),
           envir = to)
  }
}

load_model_artifact <- function(manifest_path = NULL, artifact,
                                envir = parent.frame(), legacy_path = NULL) {
  manifest_missing <- is_missing_path(manifest_path)
  legacy_missing <- is_missing_path(legacy_path)
  if (manifest_missing == legacy_missing) {
    stop(
      "Supply exactly one of manifest_path or the explicit legacy_path.",
      call. = FALSE
    )
  }

  if (!legacy_missing) {
    warning(
      "Loading an unmanifested legacy model artifact; schema, status, and gene-universe validation are unavailable.",
      call. = FALSE
    )
    inspected <- inspect_model_artifact(
      legacy_path, artifact, require_status = FALSE
    )
    copy_artifact_objects(
      inspected$environment, envir, inspected$required_objects
    )
    return(invisible(list(path = inspected$path, manifest = NULL,
                          artifact = artifact, legacy = TRUE)))
  }

  manifest <- read_model_manifest(manifest_path)
  entry <- manifest$artifacts[[artifact]]
  if (is.null(entry)) {
    stop("Manifest does not contain artifact '", artifact, "'.", call. = FALSE)
  }
  if (!is.list(entry) || !all(c("path", "required_objects", "gene_ids",
                                "fit_records") %in% names(entry))) {
    stop("Manifest entry '", artifact, "' is malformed.", call. = FALSE)
  }
  validate_character_scalar(entry$path, "Manifest artifact path")
  if (!file.exists(entry$path)) {
    stop("Manifest-selected model artifact does not exist: ", entry$path,
         call. = FALSE)
  }
  if (!identical(normalizePath(entry$path, mustWork = TRUE), entry$path)) {
    stop("Manifest artifact path is not the exact normalized path.",
         call. = FALSE)
  }
  inspected <- inspect_model_artifact(entry$path, artifact)
  if (!identical(inspected$required_objects, entry$required_objects)) {
    stop("Manifest required-object record does not match the artifact.",
         call. = FALSE)
  }
  if (!identical(inspected$gene_ids, entry$gene_ids)) {
    stop("Manifest gene universe or order does not match the artifact.",
         call. = FALSE)
  }
  if (!identical(inspected$fit_records, entry$fit_records)) {
    stop("Manifest fit-status records do not match the artifact.",
         call. = FALSE)
  }
  if (any(inspected$fit_records$optimizer != manifest$settings$optimizer)) {
    stop("Manifest optimizer does not match the artifact fits.", call. = FALSE)
  }
  if (any(inspected$fit_records$bootstrap_replicates !=
          manifest$settings$bootstrap_count)) {
    stop("Manifest bootstrap count does not match the artifact fits.",
         call. = FALSE)
  }
  if (any(inspected$fit_records$grid_size != manifest$settings$grid_size)) {
    stop("Manifest grid size does not match the artifact fits.", call. = FALSE)
  }
  effect_setting <- if (identical(artifact, "ddd")) "ddd" else "autism"
  declared_effect <- manifest$settings$max_effect_size[[effect_setting]]
  if (is.null(declared_effect) || any(abs(
    inspected$fit_records$maximum_effect_size - declared_effect
  ) > 1e-10 * max(1, declared_effect))) {
    stop("Manifest maximum effect size does not match the artifact fits.",
         call. = FALSE)
  }
  unreliable <- inspected$fit_records$fit_id[
    !is.na(inspected$fit_records$uncertainty_reliable) &
      !inspected$fit_records$uncertainty_reliable
  ]
  if (length(unreliable)) {
    warning(
      "Manifest-selected artifact has unreliable uncertainty for fits: ",
      paste(unreliable, collapse = ", "), ".", call. = FALSE
    )
  }
  copy_artifact_objects(
    inspected$environment, envir, inspected$required_objects
  )
  invisible(list(path = inspected$path, manifest = manifest,
                 artifact = artifact, legacy = FALSE))
}
