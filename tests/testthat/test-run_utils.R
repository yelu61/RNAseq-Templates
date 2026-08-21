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
    role = "canonical", retention = "full", change_note = "primary"
  )

  expect_true(file.exists(file.path(run_dir, "run_manifest.csv")))
  expect_equal(x$run_id, basename(run_dir))
  expect_equal(x$role, "canonical")
  expect_match(x$input_md5, "^[0-9a-f]{32}$")
  expect_match(x$analysis_signature, "^[0-9a-f]{32}$")
  expect_match(x$backend_signature, "^[0-9a-f]{32}$")
})

test_that("build_run_registry identifies exact duplicates but not changed configs", {
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

test_that("a canonical full run owns an exact-duplicate family", {
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
