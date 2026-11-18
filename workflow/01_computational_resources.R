## ===========================================================================
## 01_computational_resources.R - runtime, CPU time, peak memory, CPU cores and
## disk I/O per pulposeq run and per pipeline step, from the Nextflow execution
## trace (LongBench ONT cDNA, ONT DRS, PacBio; Dong ONT cDNA). Extended Data
## figure (computational resources) + the tables behind it (Supplementary Table).
##
## One trace per run: the newest execution_trace_<date>_<time>.txt, chosen by the
## date and time in its file name. Its COMPLETED and CACHED tasks are the run's
## tasks: a CACHED row keeps the metrics and the submit time of the run that
## computed it, so it is placed and timed where it ran.
## Queue time is NOT counted: the trace has submit, duration (submit -> end,
## queue included) and realtime (the task's own run time), so
##   end = submit + duration,  start = end - realtime.
## Execution time of a run = idealised critical-path (uncontended) time:
## a step takes its slowest task's realtime (all its tasks in parallel) and the
## run the longest chain of steps through its own Nextflow DAG (pipeline_dag_*.html).
## Queue time and HPC partition limits are removed; the time observed on the
## cluster (union of task intervals, queue excluded) is kept as elapsed_h.
## Per step: slowest_task_h (uncontended), realtime_h (union), realtime_sum_h.
## Per step (process name without the workflow prefix and sample tag): CPU time
## = sum of realtime x %cpu, peak memory = the largest task's peak RSS, CPU cores
## = the busiest task's %cpu / 100, data read / written = summed rchar / wchar.
## Failed attempts count towards the run's realtime and are listed: failed_attempts,
## their CPU time (cpu_h_failed) and failed_without_metrics, the attempts the trace
## records without realtime or %cpu ("-"). cpu_h_failed is NA, not 0, when no
## failed attempt has both, so a missing value is not reported as no CPU time.
## idle_h = span_h - elapsed_h: time between the first start and the last end with
## no task executing (waits between steps and, in a resumed run, the gap before
## the resume).
## Trace times are kept as the trace writes them (the cluster's clock time, read
## as UTC so that no time-zone conversion is applied on the way to the CSV and the
## workbook): first_start and last_end in Supplementary Table 2 match the trace.
##
##   a  per run              c  per-step resource profile: realtime (sqrt) vs
##   b  per pipeline step       peak RAM (log), point area = CPU cores of the
##                              busiest task, colour = pipeline stage (Fig. 1a)
##                           d  disk I/O of the main steps: data read vs written
##
## Inputs (read only): <project>/{ont_cdna,ont_drs,pb,dong}/pipeline_info/execution_trace_*.txt
##                     and pipeline_dag_*.html
## Outputs: resources/results/resources_{per_task,per_step,per_run}.csv
##          <project>/figures/"Extended data figure 1".{pdf,svg,png}
## ===========================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source("C:/Users/mbarb/OneDrive/paper-pulposeq/pulposeq-benchmarking-2026/workflow/downstream/longbench/00_paths.R")
source(file.path(repo, "downstream/longbench/helpers/plot_theme.R"))      # font_family, protocol_cols
source(file.path(repo, "downstream/longbench/helpers/figure_helpers.R"))  # save_panel, ggsave_svg
out_dir <- file.path(repo, "downstream/resources/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- c(ont_cdna = "LongBench ONT cDNA", ont_drs = "LongBench ONT DRS",
          pb = "LongBench PacBio", dong = "Dong ONT cDNA")
run_cols <- c("LongBench ONT cDNA" = protocol_cols[["ONT cDNA"]], "LongBench ONT DRS" = protocol_cols[["ONT DRS"]],
              "LongBench PacBio" = protocol_cols[["PacBio"]], "Dong ONT cDNA" = "#27408B")

## pipeline stages (the four subworkflows of Fig. 1a); a process missing here stops the script
stage_map <- list(
  "1 · QC and mapping" = c("NANOCOMP_RAW", "NANOCOMP_REST", "NANOCOMP_MINIMAP2", "RESTRANDER",
                           "MINIMAP2_ALIGN", "BAM_COVERAGE"),
  "2 · Assembly, quantification, coding potential" =
    c("BAMBU", "PREFILTER", "GFFREAD", "GFFCOMPARE", "CPC2",
      "OARFISH_INDEX_FASTA", "OARFISH_READMODE", "OARFISH_MERGE"),
  "3 · Annotation refinement, lncRNA discovery" =
    c("VALIDATE_BAMBU_OUTPUTS", "NOVEL_TRANSCRIPTS", "REF_TRANSCRIPTS", "ANNOTATION", "CURATE"),
  "4 · Visualization and reporting" = c("GENOMIC_CONTEXT", "RENDER_REPORT", "MULTIQC"))
stage_tbl  <- tibble(stage = rep(names(stage_map), lengths(stage_map)), step = unlist(stage_map, use.names = FALSE))
stage_cols <- setNames(c("#6E0B47", "#B8327F", "#E08AB9", "#F5D3E6"), names(stage_map))
step_labels <- c(NANOCOMP_MINIMAP2 = "NanoComp (BAM)", RESTRANDER = "Restrander", MINIMAP2_ALIGN = "minimap2",
                 BAM_COVERAGE = "Coverage", BAMBU = "Bambu", OARFISH_READMODE = "oarfish read mode",
                 VALIDATE_BAMBU_OUTPUTS = "Validate outputs")
io_steps <- c("RESTRANDER", "MINIMAP2_ALIGN", "BAM_COVERAGE", "BAMBU", "VALIDATE_BAMBU_OUTPUTS")

dur_s <- function(x) vapply(x, function(s) {
  if (is.na(s) || s %in% c("-", "")) return(NA_real_)
  m <- regmatches(s, gregexpr("[0-9.]+(ms|d|h|m|s)", s))[[1]]
  sum(vapply(m, function(t) { v <- as.numeric(sub("[a-z]+$", "", t)); u <- sub("^[0-9.]+", "", t)
    v * c(ms = 0.001, s = 1, m = 60, h = 3600, d = 86400)[[u]] }, 0))
}, 0, USE.NAMES = FALSE)
bytes <- function(x) vapply(x, function(s) {
  if (is.na(s) || s %in% c("-", "")) return(NA_real_)
  p <- strsplit(trimws(s), " ")[[1]]; v <- as.numeric(p[1]); if (length(p) == 1) return(v)
  v * c(B = 1, KB = 1024, MB = 1024^2, GB = 1024^3, TB = 1024^4)[[p[2]]]
}, 0, USE.NAMES = FALSE)
pct <- function(x) as.numeric(sub("%", "", ifelse(x %in% c("-", ""), NA, x)))
step_of <- function(name) sub("^.*:", "", sub(" \\(.*$", "", name))
## hours with at least one interval open (union of [start, end])
union_h <- function(s, e) {
  ok <- !is.na(s) & !is.na(e); s <- as.numeric(s[ok]); e <- as.numeric(e[ok])
  if (!length(s)) return(NA_real_)
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) if (s[i] <= ce) ce <- max(ce, e[i]) else { tot <- tot + ce - cs; cs <- s[i]; ce <- e[i] }
  (tot + ce - cs) / 3600
}
## Idealised critical-path (uncontended) time of a run: every task of
## a step runs at once, so a step takes its slowest task's realtime; steps follow
## the run's own Nextflow DAG (pipeline_dag_*.html, mermaid), and where independent
## branches run concurrently the longest one counts. Removes HPC queue and
## partition limits; measures the workflow itself.
dag_edges <- function(dir) {
  fs <- list.files(dir, "^pipeline_dag_.*\\.html$", full.names = TRUE)
  if (!length(fs)) stop("no pipeline_dag_*.html in ", dir)
  html  <- readLines(fs[order(basename(fs), decreasing = TRUE)][1], warn = FALSE)
  nodes <- unlist(regmatches(html, gregexpr('v[0-9]+\\(\\["[^"]+"\\]\\)', html)))
  proc  <- setNames(sub('^.*\\["([^"]+)"\\].*$', "\\1", nodes), sub("\\(.*$", "", nodes))
  e     <- do.call(rbind, strsplit(gsub(" ", "", unlist(regmatches(html, gregexpr("v[0-9]+ *--> *v[0-9]+", html)))), "-->"))
  succ  <- split(e[, 2], e[, 1])
  ## processes reached from a process through channel / operator nodes only
  nxt <- function(v) {
    seen <- character(0); todo <- succ[[v]]; out <- character(0)
    while (length(todo)) {
      x <- todo[1]; todo <- todo[-1]
      if (x %in% seen) next; seen <- c(seen, x)
      if (x %in% names(proc)) out <- c(out, x) else todo <- c(todo, succ[[x]])
    }
    unique(out)
  }
  bind_rows(lapply(names(proc), function(v) { r <- nxt(v)
    if (length(r)) tibble(from = proc[[v]], to = unname(proc[r])) else NULL })) %>% distinct()
}
critical_path <- function(step_h, edges) {          # step_h: named hours (slowest task per step)
  fin <- numeric(0); via <- character(0)
  f_of <- function(s) {
    if (s %in% names(fin)) return(fin[[s]])
    pr <- intersect(unique(edges$from[edges$to == s]), names(step_h))
    best <- 0
    if (length(pr)) { v <- vapply(pr, f_of, 0); best <- max(v); via[[s]] <<- pr[which.max(v)] }
    fin[[s]] <<- best + step_h[[s]]
    fin[[s]]
  }
  for (s in names(step_h)) f_of(s)
  path <- names(which.max(fin))
  while (path[1] %in% names(via)) path <- c(via[[path[1]]], path)
  list(hours = max(fin), path = path)
}

## the newest trace: by the date and time in its name (execution_trace_YYYY-MM-DD_HH-MM-SS.txt)
newest_trace <- function(dir) {
  fs <- list.files(dir, "^execution_trace_.*\\.txt$", full.names = TRUE)
  if (!length(fs)) stop("no execution trace in ", dir)
  t <- as.POSIXct(sub("^execution_trace_(.*)\\.txt$", "\\1", basename(fs)), format = "%Y-%m-%d_%H-%M-%S", tz = "UTC")
  if (anyNA(t)) stop("trace names without a date and time: ", paste(basename(fs)[is.na(t)], collapse = ", "))
  fs[which.max(t)]
}

per_task <- list(); per_run <- list()
for (p in names(runs)) {
  f  <- newest_trace(file.path(project, p, "pipeline_info"))
  tr <- read.delim(f, colClasses = "character", check.names = FALSE)
  n_cached <- sum(tr$status == "CACHED")
  if (n_cached) message(sprintf("  %s: %d of %d tasks CACHED in %s (counted with their own metrics)",
                                p, n_cached, sum(tr$status %in% c("COMPLETED", "CACHED")), basename(f)))
  ran <- tr %>% filter(status %in% c("COMPLETED", "CACHED", "FAILED")) %>%
    ## tz = "UTC": the trace's clock time, kept as written (see the header)
    mutate(submit_t = as.POSIXct(submit, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
           end      = submit_t + dur_s(duration),
           start    = end - dur_s(realtime),
           queue_h  = (dur_s(duration) - dur_s(realtime)) / 3600)
  done   <- ran %>% filter(status %in% c("COMPLETED", "CACHED")) %>%
    mutate(run = runs[[p]], trace = basename(f), step = step_of(name),
           realtime_h = dur_s(realtime) / 3600, cores = pct(`%cpu`) / 100, cpu_h = realtime_h * cores,
           peak_rss_gb = bytes(peak_rss) / 1024^3, read_gb = bytes(rchar) / 1024^3, written_gb = bytes(wchar) / 1024^3)
  failed <- ran %>% filter(status == "FAILED") %>%
    mutate(cpu_h = dur_s(realtime) / 3600 * pct(`%cpu`) / 100)   # NA when the trace has "-"
  if (any(is.na(failed$cpu_h)))
    message(sprintf("  %s: %d of %d failed attempts have no realtime or %%cpu in the trace",
                    p, sum(is.na(failed$cpu_h)), nrow(failed)))
  per_task[[p]] <- done %>% select(run, trace, step, name, hash, start, end, realtime_h, queue_h, cores, cpu_h,
                                   peak_rss_gb, read_gb, written_gb)
  slowest <- done %>% group_by(step) %>% summarise(h = max(realtime_h, na.rm = TRUE), .groups = "drop")
  cp <- critical_path(setNames(slowest$h, slowest$step), dag_edges(file.path(project, p, "pipeline_info")))
  top <- done %>% slice_max(peak_rss_gb, n = 1, with_ties = FALSE)
  busy <- done %>% slice_max(cores, n = 1, with_ties = FALSE)
  per_run[[p]] <- tibble(
    run = runs[[p]], trace = basename(f), samples = sum(done$step == "MINIMAP2_ALIGN"),
    final_tasks = nrow(done), cached_tasks = n_cached, failed_attempts = nrow(failed),
    failed_without_metrics = sum(is.na(failed$cpu_h)),
    critical_path_h = cp$hours,  # idealised, uncontended
    critical_path = paste(cp$path, collapse = " > "),
    elapsed_h = union_h(ran$start, ran$end),                                   # observed on the HPC, queue excluded
    span_h = as.numeric(difftime(max(ran$end, na.rm = TRUE), min(ran$start, na.rm = TRUE), units = "hours")),
    idle_h = span_h - elapsed_h,                                               # no task executing
    elapsed_h_with_queue = as.numeric(difftime(max(ran$end, na.rm = TRUE), min(ran$submit_t, na.rm = TRUE), units = "hours")),
    cpu_h = sum(done$cpu_h, na.rm = TRUE),
    ## NA, not 0, when no failed attempt has realtime and %cpu in the trace
    cpu_h_failed = if (nrow(failed) && all(is.na(failed$cpu_h))) NA_real_ else sum(failed$cpu_h, na.rm = TRUE),
    peak_rss_gb = top$peak_rss_gb, peak_step = top$step, peak_cores = busy$cores, peak_cores_step = busy$step,
    read_tb = sum(done$read_gb, na.rm = TRUE) / 1024, written_tb = sum(done$written_gb, na.rm = TRUE) / 1024,
    queue_h_summed = sum(done$queue_h, na.rm = TRUE))
}
per_task <- bind_rows(per_task); per_run <- bind_rows(per_run)
per_step <- per_task %>% group_by(run, step) %>%
  summarise(tasks = n(), first_start = min(start), last_end = max(end),
            slowest_task_h = max(realtime_h, na.rm = TRUE),  # uncontended step time
            ## the sum, before realtime_h is redefined as the union of task intervals
            realtime_sum_h = sum(realtime_h), realtime_h = union_h(start, end), cpu_h = sum(cpu_h, na.rm = TRUE),
            max_cores = max(cores, na.rm = TRUE), peak_rss_gb = max(peak_rss_gb, na.rm = TRUE),
            read_gb = sum(read_gb, na.rm = TRUE), written_gb = sum(written_gb, na.rm = TRUE), .groups = "drop") %>%
  left_join(stage_tbl, by = "step")
unmapped <- unique(per_step$step[is.na(per_step$stage)])
if (length(unmapped)) stop("processes without a stage in stage_map: ", paste(unmapped, collapse = ", "))
write_csv(per_task, file.path(out_dir, "resources_per_task.csv"))
write_csv(per_step, file.path(out_dir, "resources_per_step.csv"))
write_csv(per_run,  file.path(out_dir, "resources_per_run.csv"))
print(per_run %>% select(run, trace, samples, final_tasks, cached_tasks, failed_attempts, elapsed_h, span_h,
                         cpu_h, peak_rss_gb, peak_step, peak_cores, written_tb), width = 220)


## ---------------------------------------------------------------- figure ---
## Laid out at print size for the A4 text block under the caption (about
## 155 x 185 mm): every size below is the printed size. Text is 6-7 pt,
## panel letters 9 pt. Change fig_w / fig_h if the page margins change.
##
##   a | b          a  per run (bars, value printed)        b  per step (log)
##   --+---------   c  per-step profile, one facet per run (2 x 2)
##   c     | d      d  disk I/O of the main steps (colours as in a and b)
##         | key    key: pipeline stage and CPU-core legends of c
fig_w   <- 200   # mm, figure width (wide enough for one-line strip titles in b)
fig_h   <- 185   # mm, height left under the caption
top_h   <- 72    # mm, row of a and b
a_w     <- 40    # mm, column of a
right_w <- 58    # mm, column of d and the legend of c
d_h     <- 60    # mm, frame of d (square data area + axes; its shape legend sits inside)
gap     <- 3     # mm between rows and columns
tag_mm  <- 3.5   # mm left margin of every panel: room for its letter
fs_txt  <- 6     # pt: tick labels, legend entries, point labels
fs_ttl  <- 7     # pt: axis titles, facet strips, legend titles
fs_tag  <- 9     # pt: panel letters

theme_set(theme_bw(base_family = font_family, base_size = fs_ttl) +
  theme(text = element_text(face = "plain"), axis.text = element_text(size = fs_txt, colour = "grey15"),
        axis.title = element_text(size = fs_ttl),
        strip.text = element_text(size = fs_ttl, margin = margin(0.6, 0.6, 0.6, 0.6, "mm")),
        strip.background = element_rect(fill = "grey94", linewidth = 0.25), panel.border = element_rect(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25), axis.ticks.length = unit(0.6, "mm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        legend.text = element_text(size = fs_txt), legend.title = element_text(size = fs_ttl),
        legend.key.size = unit(2.6, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
        plot.margin = margin(0.5, 0.5, 0.5, tag_mm, "mm")))
run_lv    <- unname(runs)
run_short <- c("LongBench ONT cDNA" = "ONT cDNA", "LongBench ONT DRS" = "ONT DRS",
               "LongBench PacBio" = "PacBio", "Dong ONT cDNA" = "ONT cDNA (Dong)")

## a: per run, horizontal bars with the value printed (no x axis needed)
ra <- per_run %>% transmute(run = factor(run, levels = run_lv),
                            `Execution time (h)` = critical_path_h, `CPU time (h)` = cpu_h,
                            `Peak memory (GB)` = peak_rss_gb, `Data written (TB)` = written_tb) %>%
  pivot_longer(-run, names_to = "metric") %>%
  mutate(metric = factor(metric, levels = unique(metric)),
         lab = ifelse(value >= 100, sprintf("%.0f", value), ifelse(value >= 10, sprintf("%.1f", value), sprintf("%.2f", value))))
pA <- ggplot(ra, aes(value, fct_rev(run), fill = run)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = lab), hjust = -0.15, size = fs_txt / .pt, family = font_family) +
  facet_wrap(~ metric, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = run_cols, guide = "none") +
  scale_y_discrete(labels = run_short) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), panel.spacing.y = unit(1, "mm"))

## b: per step (dots, log scale; steps in pipeline order = median relative first start)
step_ord <- per_step %>% group_by(run) %>% mutate(rel = as.numeric(difftime(first_start, min(first_start), units = "secs"))) %>%
  group_by(step) %>% summarise(rel = median(rel), .groups = "drop") %>% arrange(rel)
rb <- per_step %>% transmute(run = factor(run, levels = run_lv), step = factor(step, levels = rev(step_ord$step)),
                             `Slowest task (h)` = slowest_task_h, `CPU time (h)` = cpu_h,
                             `Peak memory (GB)` = peak_rss_gb, `Data written (GB)` = written_gb) %>%
  pivot_longer(-c(run, step), names_to = "metric") %>%
  mutate(metric = factor(metric, levels = unique(metric)), value = pmax(value, 1e-3))
pB <- ggplot(rb, aes(value, step, colour = run)) +
  geom_point(size = 0.8, alpha = 0.85, position = position_dodge(width = 0.6)) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  ## one break per decade, every other decade when a facet spans more than four
  scale_x_log10(breaks = function(l) { e <- seq(ceiling(log10(l[1])), floor(log10(l[2])))
                                       if (length(e) > 4) e <- e[(max(e) - e) %% 2 == 0]; 10^e },
                labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE, trim = TRUE),
                guide = guide_axis(check.overlap = TRUE)) +   # drop tick labels that would collide
  scale_colour_manual(values = run_cols, labels = run_short, name = NULL,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "bottom", legend.key.spacing.x = unit(1.5, "mm"),
        panel.spacing.x = unit(4, "mm"))

## c: per-step resource profile, one facet per run (2 x 2)
size_range <- c(0.4, 3.6)
strip_lab <- setNames(sprintf("%s (n = %d)", per_run$run, per_run$samples), per_run$run)
prof <- per_step %>% mutate(run = factor(run, levels = run_lv), stage = factor(stage, levels = names(stage_map)),
                            label = unname(step_labels[step]),
                            pt_size = size_range[1] + sqrt(scales::rescale(max_cores)) * diff(size_range))
x_max <- max(24, ceiling(max(prof$slowest_task_h, na.rm = TRUE)))
y_lim <- c(min(5e-3, min(prof$peak_rss_gb, na.rm = TRUE)), max(5000, max(prof$peak_rss_gb, na.rm = TRUE) * 5))
pC_full <- ggplot(prof, aes(slowest_task_h, peak_rss_gb)) +
  geom_point(aes(size = max_cores, fill = stage), shape = 21, colour = "black", stroke = 0.2, alpha = 0.9) +
  ## every point is passed (unlabelled ones as ""), so the labels also keep clear of them
  ggrepel::geom_text_repel(data = mutate(prof, label = coalesce(label, "")), aes(label = label, point.size = pt_size),
                           size = fs_txt / .pt, family = font_family, min.segment.length = 0, segment.size = 0.15,
                           segment.colour = "grey30", box.padding = 0.3, point.padding = 0.25, force = 4,
                           max.overlaps = Inf, max.time = 3, max.iter = 1e5, seed = 1) +
  facet_wrap(~ run, ncol = 2, labeller = as_labeller(strip_lab)) +
  scale_x_continuous(transform = "sqrt", breaks = c(1 / 60, 1, 4, 8, 12, 24),
                     labels = c("1m", "1h", "4h", "8h", "12h", "24h"), limits = c(0, x_max)) +
  scale_y_log10(breaks = 10^(-2:3), labels = c("0.01", "0.1", "1", "10", "100", "1,000"), limits = y_lim) +
  scale_size(name = "CPU cores of the busiest task", range = size_range, breaks = c(1, 5, 10, 20),
             guide = guide_legend(order = 2, nrow = 1, direction = "horizontal", title.position = "top")) +
  scale_fill_manual(values = stage_cols, name = "Pipeline stage", drop = FALSE,
                    guide = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2.2))) +
  labs(x = "Slowest task (square-root scale)", y = "Peak RAM, largest task (GB, log)") +
  ## the legends are drawn on their own, as the key under d (legend_grob)
  theme(panel.grid.major = element_blank(), panel.spacing = unit(1.5, "mm"),
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.justification = c(0, 1), legend.key.spacing.y = unit(0.6, "mm"),
        legend.key.spacing.x = unit(1.5, "mm"), legend.title = element_text(size = fs_txt, margin = margin(b = 1, unit = "mm")),
        legend.spacing.y = unit(3.5, "mm"))
pC <- pC_full + theme(legend.position = "none")

## d: disk I/O of the main steps (equal log axes; dashed: written = read; colours as in a/b)
io <- per_step %>% filter(step %in% io_steps) %>%
  mutate(run = factor(run, levels = run_lv), step_lab = factor(step_labels[step], levels = step_labels[io_steps]))
io_lim <- c(0.1, 3000)
if (any(unlist(io[c("read_gb", "written_gb")]) < io_lim[1] | unlist(io[c("read_gb", "written_gb")]) > io_lim[2], na.rm = TRUE))
  message("  ! some disk I/O values fall outside ", io_lim[1], "-", io_lim[2], " GB (clipped): widen io_lim")
pD <- ggplot(io, aes(read_gb, written_gb)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey55", linewidth = 0.25) +
  annotate("text", x = 4, y = 11, label = "written = read", angle = 45, colour = "grey40", size = fs_txt / .pt,
           family = font_family) +
  annotate("text", x = 10, y = 0.22, label = "reads far more\nthan it writes", lineheight = 0.9,
           colour = "grey40", size = fs_txt / .pt, family = font_family) +
  geom_point(aes(fill = run, shape = step_lab), size = 1.6, colour = "black", stroke = 0.25) +
  scale_x_log10(breaks = 10^(-1:3), labels = c("0.1", "1", "10", "100", "1,000")) +
  scale_y_log10(breaks = 10^(-1:3), labels = c("0.1", "1", "10", "100", "1,000")) +
  coord_fixed(xlim = io_lim, ylim = io_lim, expand = FALSE) +
  scale_fill_manual(values = run_cols, guide = "none") +
  scale_shape_manual(values = setNames(c(21, 22, 24, 23, 25), step_labels[io_steps]), name = "Step",
                     guide = guide_legend(ncol = 1, override.aes = list(fill = "grey85", size = 1.8))) +
  labs(x = "Data read (GB, log)", y = "Data written (GB, log)") +
  ## step legend inside the empty upper-left corner (above the diagonal)
  theme(legend.position = "inside", legend.position.inside = c(0.03, 0.97), legend.justification.inside = c(0, 1),
        legend.background = element_blank(), legend.key = element_blank(), legend.key.spacing.y = unit(0.3, "mm"),
        legend.title = element_text(size = fs_txt), panel.grid.major = element_blank())

## Two helpers drawn at render time, so text is measured on the output device
## (showtext / svglite), exactly as patchwork does for a ggplot:
##   square_grob(): a coord_fixed plot as the largest square that fits its frame,
##                  pinned to the top-left (coord_fixed would otherwise centre it)
##   legend_grob(): only the legends of a plot, pinned to the top-left
square_grob <- function(p) grid::gTree(plot = p, cl = "square_plot")
registerS3method("makeContent", "square_plot", function(x) {
  g  <- ggplotGrob(x$plot)
  pn <- g$layout[grepl("^panel", g$layout$name), ]
  mm_h <- function(u) grid::convertHeight(u, "mm", valueOnly = TRUE)
  mm_w <- function(u) grid::convertWidth(u, "mm", valueOnly = TRUE)
  side <- min(mm_h(unit(1, "npc")) - mm_h(sum(g$heights[-pn$t])),
              mm_w(unit(1, "npc")) - mm_w(sum(g$widths[-pn$l])))
  g$widths[pn$l]  <- unit(side, "mm")
  g$heights[pn$t] <- unit(side, "mm")
  g$respect <- FALSE
  g <- gtable::gtable_add_cols(g, unit(1, "null"))   # free space to the right
  g <- gtable::gtable_add_rows(g, unit(1, "null"))   # free space below
  grid::setChildren(x, grid::gList(g))
}, envir = asNamespace("grid"))
legend_grob <- function(p) grid::gTree(plot = p, cl = "legend_only")
registerS3method("makeContent", "legend_only", function(x) {
  g   <- ggplotGrob(x$plot)
  i   <- grep("^guide-box", g$layout$name)
  leg <- g$grobs[i][!vapply(g$grobs[i], inherits, TRUE, "zeroGrob")][[1]]
  leg$vp <- grid::viewport(x = unit(tag_mm, "mm"), y = 1, just = c("left", "top"),
                           width = grid::grobWidth(leg), height = grid::grobHeight(leg))
  grid::setChildren(x, grid::gList(leg))
}, envir = asNamespace("grid"))

## grid: [a | gap | middle | gap | right] x [top | gap | d | key]
design <- "A#BBB
           #####
           CCC#D
           CCC#E"
p <- wrap_elements(full = pA) + wrap_elements(full = pB) + wrap_elements(full = pC) +
  wrap_elements(full = square_grob(pD)) + wrap_elements(full = legend_grob(pC_full)) +
  plot_layout(design = design,
              widths  = unit(c(a_w, gap, 1, gap, right_w), c("mm", "mm", "null", "mm", "mm")),
              heights = unit(c(top_h, gap, d_h, 1), c("mm", "mm", "mm", "null"))) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", ""))) &   # the key (E) has no letter
  theme(plot.tag = element_text(size = fs_tag, face = "plain", family = font_family),
        plot.margin = margin(0.5, 0.5, 0.5, tag_mm, "mm"))
save_panel(p, "Extended data figure 1", c(fig_w, fig_h), file.path(project, "figures"))
message("Done: ", out_dir)
