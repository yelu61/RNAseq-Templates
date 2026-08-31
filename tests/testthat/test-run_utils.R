test_that("write_run_manifest records checksums and lifecycle metadata", {
  run_dir <- tempfile("run-")
  dir.create(run_dir)
  input <- file.path(run_dir, "counts.tsv")
  config <- file.path(run_dir, "config.R")
  writeLines("gene\tS1\nA\t1", input)
  writeLines("SPECIES <- 'mouse'", config)

  x <- write_run_manifest(
    run_dir, config, input,
    analysis_config = list(SPECIES = "mouse", MIN_COUNT = 10),
    role = "canonical", retention = "full", change_note = "primary",
    input_files = character(0), backend_files = config
  )

  expect_true(file.exists(file.path(run_dir, "run_manifest.csv")))
  expect_equal(x$run_id, basename(run_dir))
  expect_equal(x$role, "canonical")
  expect_match(x$input_md5, "^[0-9a-f]{32}$")
  expect_match(x$analysis_signature, "^[0-9a-f]{32}$")
  expect_match(x$backend_signature, "^[0-9a-f]{32}$")
  expect_equal(x$manifest_schema_version, 2L)
  expect_true(x$inputs_complete)
  expect_true(x$duplicate_eligible)
  runtime <- read.csv(file.path(run_dir, "run_runtime.csv"))
  expect_equal(runtime$version[runtime$component == "R"], R.version$version.string)
  expect_match(x$runtime_signature, "^[0-9a-f]{32}$")
})

test_that("initialize_run_lifecycle supplies defaults and validates overrides", {
  env <- new.env(parent = emptyenv())
  expect_invisible(initialize_run_lifecycle(env))
  expect_equal(env$RUN_ROLE, "candidate")
  expect_equal(env$RUN_RETENTION, "full")

  env$RUN_ROLE <- "latest"
  expect_error(initialize_run_lifecycle(env), "RUN_ROLE")
})

test_that("build_run_registry groups complete matching provenance but not changed configs", {
  runs <- tempfile("runs-")
  dir.create(runs)
  input <- tempfile(fileext = ".tsv")
  writeLines("gene\tS1\nA\t1", input)
  for (i in 1:3) {
    run <- file.path(runs, paste0("run", i))
    dir.create(run)
    config <- file.path(run, "config.R")
    writeLines("SPECIES <- 'mouse'", config)
    cfg <- list(SPECIES = "mouse", MIN_COUNT = if (i == 3) 20 else 10)
    write_run_manifest(run, config, input, cfg,
                       input_files = character(0), backend_files = config,
                       completed_at = as.POSIXct(sprintf("2026-01-0%d 00:00:00", i), tz = "UTC"))
  }

  registry <- build_run_registry(runs)
  expect_equal(registry$duplicate_of, c("", "run1", ""))
  expect_true(file.exists(file.path(runs, "RUN_REGISTRY.csv")))
})

test_that("run manifest rejects unsupported lifecycle values", {
  run <- tempfile("run-")
  dir.create(run)
  input <- file.path(run, "x")
  writeLines("x", input)
  expect_error(write_run_manifest(run, input, input, list(a = 1), role = "latest"),
               "Unsupported run role")
  expect_error(write_run_manifest(run, input, input, list(a = 1), retention = "delete"),
               "Unsupported run retention")
})

test_that("analysis signatures change with gene sets but not lifecycle metadata", {
  base <- list(SPECIES = "mouse", custom_gene_sets = list(IFN = c("Stat1", "Irf7")))
  changed <- list(SPECIES = "mouse", custom_gene_sets = list(IFN = c("Stat1", "Irf7", "Isg15")))
  expect_false(identical(.run_md5_object(base), .run_md5_object(changed)))

  with_lifecycle <- c(base, list(RUN_ROLE = "sensitivity"))
  expect_identical(.run_md5_object(base), .run_md5_object(with_lifecycle[names(base)]))
})

test_that("formula signatures do not serialize their environments", {
  env1 <- new.env(parent = emptyenv())
  env2 <- new.env(parent = emptyenv())
  f1 <- as.formula("~ batch + condition", env = env1)
  f2 <- as.formula("~ batch + condition", env = env2)
  env1$large_unrelated_object <- raw(10000)
  expect_identical(.run_md5_object(list(design = f1)),
                   .run_md5_object(list(design = f2)))
})

test_that("a canonical full run owns a matching-provenance family", {
  runs <- tempfile("runs-")
  dir.create(runs)
  input <- tempfile()
  writeLines("counts", input)
  roles <- c("candidate", "canonical", "repro_check")
  retention <- c("full", "full", "metadata_only")
  for (i in seq_along(roles)) {
    run <- file.path(runs, paste0("run", i))
    dir.create(run)
    config <- file.path(run, "config.R")
    writeLines("same", config)
    write_run_manifest(run, config, input, list(a = 1), role = roles[[i]],
                       input_files = character(0), backend_files = config,
                       retention = retention[[i]], completed_at = as.POSIXct(i, origin = "2026-01-01", tz = "UTC"))
  }
  registry <- build_run_registry(runs)
  expect_equal(registry$duplicate_of[registry$run_id == "run2"], "")
  expect_true(all(registry$duplicate_of[registry$run_id != "run2"] == "run2"))
})

test_that("backend signatures change when analysis code changes", {
  code <- tempfile(fileext = ".R")
  writeLines("x <- 1", code)
  first <- .run_backend_signature(code, revision = "abc")
  writeLines("x <- 2", code)
  second <- .run_backend_signature(code, revision = "abc")
  expect_false(identical(first, second))
  unlink(code)
})

test_that("run registry accepts legacy and current manifest schemas", {
  runs <- tempfile("runs-")
  dir.create(runs)
  input <- tempfile()
  writeLines("counts", input)
  for (id in c("legacy", "current")) {
    run <- file.path(runs, id)
    dir.create(run)
    config <- file.path(run, "config.R")
    writeLines("same", config)
    write_run_manifest(run, config, input, list(a = 1))
  }
  legacy_file <- file.path(runs, "legacy", "run_manifest.csv")
  legacy <- read.csv(legacy_file, stringsAsFactors = FALSE, check.names = FALSE)
  legacy <- legacy[, !colnames(legacy) %in% c("backend_revision", "backend_dirty", "backend_signature"), drop = FALSE]
  write.csv(legacy, legacy_file, row.names = FALSE)

  registry <- build_run_registry(runs)
  expect_equal(nrow(registry), 2)
  expect_true("backend_signature" %in% colnames(registry))
  expect_true(is.na(registry$backend_signature[registry$run_id == "legacy"]))
  expect_true(all(registry$duplicate_of == ""))
})

test_that("a primary-only manifest cannot hide changed auxiliary inputs", {
  base <- tempfile("auxiliary-inputs-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  input <- file.path(base, "counts.csv")
  metadata <- file.path(base, "metadata.csv")
  config <- file.path(base, "config.R")
  writeLines("counts", input)
  writeLines("same config", config)
  runs <- file.path(base, "runs")
  dir.create(runs)
  for (id in c("before", "after")) {
    run <- file.path(runs, id)
    dir.create(run)
    writeLines(id, metadata)
    write_run_manifest(run, config, input, list(META_FILE = metadata))
  }
  expect_true(all(build_run_registry(runs)$duplicate_of == ""))
})

test_that("unknown online input checksums never establish duplicates", {
  runs <- tempfile("online-runs-")
  dir.create(runs)
  on.exit(unlink(runs, recursive = TRUE), add = TRUE)
  for (id in c("online1", "online2")) {
    run <- file.path(runs, id)
    dir.create(run)
    config <- file.path(run, "config.R")
    writeLines("same config", config)
    write_run_manifest(run, config, NA_character_, list(GEO_ACCESSION = "same"),
                       input_files = character(0), backend_files = config)
  }
  expect_true(all(build_run_registry(runs)$duplicate_of == ""))
})

test_that("auxiliary metadata and reference contents participate in input signatures", {
  base <- tempfile("complete-inputs-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  paths <- file.path(base, c("counts.csv", "metadata.csv", "reference.gmt", "config.R", "backend.R"))
  for (path in paths) writeLines("initial", path)
  runs <- file.path(base, "runs")
  dir.create(runs)
  ids <- c("baseline", "repeat", "metadata_changed", "reference_changed")
  for (i in seq_along(ids)) {
    run <- file.path(runs, ids[[i]])
    dir.create(run)
    if (i == 3) writeLines("changed groups", paths[[2]])
    if (i == 4) writeLines("changed gene set", paths[[3]])
    auxiliary <- c(metadata = paths[[2]], reference = paths[[3]])
    if (i == 2) auxiliary <- rev(auxiliary)
    write_run_manifest(run, paths[[4]], paths[[1]], list(META_FILE = paths[[2]]),
                       input_files = auxiliary, backend_files = paths[[5]],
                       completed_at = as.POSIXct(i, origin = "2026-01-01", tz = "UTC"))
  }
  registry <- build_run_registry(runs)
  expect_equal(registry$duplicate_of, c("", "baseline", "", ""))
  expect_true(all(registry$duplicate_eligible))
  expect_length(unique(registry$input_signature), 3L)
  expect_length(unique(registry$analysis_signature), 1L)
  inventory <- read.csv(file.path(runs, "baseline", "run_inputs.csv"))
  expect_setequal(inventory$role, c("primary", "metadata", "reference"))
  expect_true(all(inventory$status == "present"))
  expect_true(all(grepl("^[0-9a-f]{32}$", inventory$md5)))
})

test_that("missing, unknown and explicitly incomplete references disable duplicate matching", {
  base <- tempfile("incomplete-inputs-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  input <- file.path(base, "counts.csv")
  writeLines("counts", input)
  for (kind in c("missing", "unknown", "incomplete", "unknown_backend")) {
    run <- file.path(base, kind)
    dir.create(run)
    auxiliary <- switch(kind, missing = c(reference = file.path(base, "absent.gmt")),
                        unknown = c(reference = NA_character_), character(0))
    manifest <- write_run_manifest(
      run, input, input, list(a = 1), input_files = auxiliary,
      inputs_complete = kind != "incomplete",
      backend_files = if (kind == "unknown_backend") character(0) else input
    )
    expect_false(manifest$duplicate_eligible, info = kind)
    if (kind != "unknown_backend") expect_false(manifest$inputs_complete, info = kind)
  }
  registry <- build_run_registry(base)
  expect_true(all(!registry$duplicate_eligible))
  expect_true(all(registry$duplicate_of == ""))
})

test_that("legacy manifests stay readable without inheriting schema-v2 duplicate eligibility", {
  runs <- tempfile("legacy-runs-")
  dir.create(runs)
  on.exit(unlink(runs, recursive = TRUE), add = TRUE)
  input <- file.path(runs, "counts.csv")
  writeLines("same input and fixture code", input)
  for (id in c("legacy1", "legacy2", "current1", "current2")) {
    run <- file.path(runs, id)
    dir.create(run)
    manifest <- write_run_manifest(run, input, input, list(a = 1),
                                   input_files = character(0), backend_files = input)
    if (startsWith(id, "legacy")) {
      manifest$manifest_schema_version <- NULL
      manifest$inputs_complete <- NULL
      manifest$input_signature <- NULL
      manifest$runtime_signature <- NULL
      manifest$duplicate_eligible <- NULL
      write.csv(manifest, file.path(run, "run_manifest.csv"), row.names = FALSE)
    }
  }
  registry <- build_run_registry(runs)
  legacy <- startsWith(registry$run_id, "legacy")
  expect_true(all(!registry$duplicate_eligible[legacy]))
  expect_true(all(registry$duplicate_of[legacy] == ""))
  expect_true(all(registry$duplicate_eligible[!legacy]))
  expect_equal(sum(nzchar(registry$duplicate_of[!legacy])), 1L)
})

test_that("different runtime fingerprints cannot be grouped", {
  runs <- tempfile("runtime-runs-")
  dir.create(runs)
  on.exit(unlink(runs, recursive = TRUE), add = TRUE)
  input <- file.path(runs, "input")
  writeLines("same", input)
  for (id in c("runtime1", "runtime2")) {
    run <- file.path(runs, id)
    dir.create(run)
    manifest <- write_run_manifest(run, input, input, list(a = 1),
                                   input_files = character(0), backend_files = input)
    if (id == "runtime2") {
      manifest$runtime_signature <- .run_md5_object(list(R = "different version"))
      write.csv(manifest, file.path(run, "run_manifest.csv"), row.names = FALSE)
    }
  }
  expect_true(all(build_run_registry(runs)$duplicate_of == ""))
})

test_that("completed run provenance cannot be overwritten", {
  run <- tempfile("protected-manifest-")
  dir.create(run)
  on.exit(unlink(run, recursive = TRUE), add = TRUE)
  input <- file.path(run, "counts.csv")
  writeLines("counts", input)
  write_run_manifest(run, input, input, list(a = 1))
  files <- file.path(run, c("run_manifest.csv", "run_inputs.csv", "run_runtime.csv"))
  before <- tools::md5sum(files)
  expect_error(write_run_manifest(run, input, input, list(a = 2)), "Refusing to overwrite")
  expect_identical(tools::md5sum(files), before)
})

test_that("fresh run checks allow input data, configuration, scripts and an open log", {
  run <- tempfile("fresh-run-")
  expect_invisible(assert_fresh_run_dir(run))
  expect_false(dir.exists(run))
  dir.create(run)
  on.exit(unlink(run, recursive = TRUE), add = TRUE)
  for (dir in c("0-Data", "config", "scripts")) dir.create(file.path(run, dir))
  for (file in c("config.R", "run_analysis.R", "run.log", "0-Data/counts.csv")) {
    writeLines("untouched", file.path(run, file))
  }
  files <- list.files(run, recursive = TRUE, full.names = TRUE)
  before <- tools::md5sum(files)
  expect_invisible(assert_fresh_run_dir(run))
  expect_identical(tools::md5sum(files), before)
})

test_that("fresh run checks reject stale optional outputs and config snapshots without mutation", {
  base <- tempfile("stale-runs-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  artifacts <- c("run_manifest.csv", "0-Config/analysis_config_used.R", "2-GSEA/old.csv",
                 "4-TME/old.csv", "5-WGCNA/old.rds", "sessionInfo.txt", "RNAseq_report.html")
  for (i in seq_along(artifacts)) {
    run <- file.path(base, i)
    file <- file.path(run, artifacts[[i]])
    dir.create(dirname(file), recursive = TRUE)
    writeLines("old output must survive", file)
    before <- tools::md5sum(file)
    expect_error(assert_fresh_run_dir(run), "Refusing to reuse", info = artifacts[[i]])
    expect_identical(tools::md5sum(file), before)
  }
})

test_that("custom output directories receive the same protection as numbered outputs", {
  base <- tempfile("custom-output-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  output <- file.path(base, "coexpression")
  expect_invisible(assert_fresh_run_dir(base, output_dirs = output))
  dir.create(output)
  expect_error(assert_fresh_run_dir(base, output_dirs = output), "Refusing to reuse")
  result <- file.path(output, "Hub_genes_all_modules.csv")
  writeLines("existing scientific result", result)
  before <- tools::md5sum(result)
  expect_error(assert_fresh_run_dir(base, output_dirs = output), "Refusing to reuse")
  expect_identical(tools::md5sum(result), before)
})

test_that("each auxiliary input content is part of the signature", {
  base <- tempfile("all-input-roles-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  primary <- file.path(base, "counts.csv")
  writeLines("counts", primary)
  roles <- c("metadata", "traits", "clinical", "tpm", "gene_sets", "annotation", "reference", "cache")
  auxiliary <- setNames(file.path(base, roles), roles)
  for (path in auxiliary) writeLines("initial", path)
  write_manifest <- function(id) {
    run <- file.path(base, id)
    dir.create(run)
    write_run_manifest(run, primary, primary, list(a = 1),
                       input_files = auxiliary, backend_files = primary)
  }
  baseline <- write_manifest("baseline")
  for (role in roles) {
    writeLines("changed", auxiliary[[role]])
    changed <- write_manifest(paste0(role, "_changed"))
    expect_true(changed$duplicate_eligible, info = role)
    expect_false(identical(changed$input_signature, baseline$input_signature), info = role)
    expect_identical(changed$analysis_signature, baseline$analysis_signature)
    writeLines("initial", auxiliary[[role]])
  }
  expect_true(all(build_run_registry(base)$duplicate_of == ""))
})

test_that("unknown or malformed checksum fields cannot establish matching provenance", {
  base <- tempfile("invalid-signatures-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  input <- file.path(base, "input")
  writeLines("same", input)
  columns <- c("input_md5", "input_signature", "config_md5", "analysis_signature",
               "backend_signature", "runtime_signature")
  for (column in columns) {
    for (value in c(NA_character_, "", "not-a-checksum")) {
      id <- paste(column, ifelse(is.na(value), "NA", ifelse(nzchar(value), "malformed", "empty")), sep = "_")
      run <- file.path(base, id)
      dir.create(run)
      manifest <- write_run_manifest(run, input, input, list(a = 1),
                                     input_files = character(0), backend_files = input)
      manifest[[column]] <- value
      write.csv(manifest, file.path(run, "run_manifest.csv"), row.names = FALSE, na = "")
    }
  }
  registry <- build_run_registry(base)
  expect_true(all(!registry$duplicate_eligible))
  expect_true(all(registry$duplicate_of == ""))
})

test_that("primary input names cannot change the reserved inventory role", {
  input <- tempfile("named-primary-")
  on.exit(unlink(input), add = TRUE)
  writeLines("counts", input)
  inventory <- .run_input_inventory(c(counts = input), character(0), dirname(input))
  expect_identical(inventory$role, "primary")
})
