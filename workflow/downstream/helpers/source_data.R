## Source Data helpers: the panel scripts save the data behind every
## ggplot they draw (the plot's data and each layer's own data), so that the
## figure data can be shared, one table per panel, without recomputing anything.
plot_data <- function(p) {
  pl <- if (inherits(p, "patchwork")) c(p$patches$plots, list(p)) else list(p)
  out <- list()
  for (q in pl) {
    if (!inherits(q, "ggplot")) next
    cand <- c(list(q$data), lapply(q$layers, function(l) l$data))
    for (d in cand) if (is.data.frame(d) && nrow(d)) {
      d <- as.data.frame(d)
      for (k in names(d)) if (is.list(d[[k]]) || is.factor(d[[k]])) d[[k]] <- vapply(d[[k]], function(v) paste(as.character(v), collapse = ";"), "")
      out[[length(out) + 1]] <- d
    }
  }
  out[!duplicated(lapply(out, function(d) c(dim(d), names(d), utils::head(unlist(d[1, ]), 20))))]
}
save_source_data <- function(plots, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(lapply(plots, function(x) plot_data(if (is.list(x) && !inherits(x, "ggplot") && !is.null(x$p)) x$p else x)), file)
  message("  source data: ", file)
}
