test_that("CLI entry points restore space-mangled Rscript paths", {
  root <- if (exists("repo_root", inherits = TRUE)) {
    get("repo_root", inherits = TRUE)
  } else {
    rprojroot::find_root(rprojroot::is_git_root, path = getwd())
  }
  files <- c(
    file.path("tools", c("build_run_registry.R", "render_report.R")),
    unlist(lapply(c("General", "Limma_Voom", "WGCNA", "TME", "TimeCourse", "TCGA_GEO"),
                  function(topic) file.path("templates", topic,
                                             c("run_analysis.R", "visualize_results.R"))))
  )
  for (file in files) {
    lines <- readLines(file.path(root, file), warn = FALSE)
    expect_true(any(grepl('gsub("~\\\\+~", " ", file_arg)', lines, fixed = TRUE)),
                info = file)
  }
})
