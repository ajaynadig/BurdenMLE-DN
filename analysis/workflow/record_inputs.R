root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"), mustWork = TRUE)
files <- sort(list.files(file.path(root, "inputs"), recursive = TRUE, full.names = TRUE))
files <- files[!file.info(files)$isdir]
manifest <- data.frame(path = substring(files, nchar(root) + 2L),
                       bytes = file.info(files)$size,
                       sha256 = vapply(files, digest::digest, character(1),
                                       algo = "sha256", file = TRUE))
dir.create(file.path(root, "provenance"), showWarnings = FALSE)
destination <- file.path(root, "provenance/input_manifest.tsv")
if (file.exists(destination)) {
  old <- read.delim(destination, check.names = FALSE)
  stopifnot(isTRUE(all.equal(old, manifest, check.attributes = FALSE)))
} else write.table(manifest, destination, sep = "\t", row.names = FALSE, quote = FALSE)
cat(nrow(manifest), "input files hashed/verified.\n")
