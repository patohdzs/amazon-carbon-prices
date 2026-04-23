#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
  cat(
    "Describe how two .Rdata files differ, including distribution shifts.\n\n",
    "Usage:\n",
    "  Rscript compare_combined_rdata.R [file_a] [file_b] [options]\n\n",
    "Defaults:\n",
    "  file_a = data/calibration/carbon_capture_curves/combined_df.Rdata\n",
    "  file_b = data/backup/calibration/carbon_capture_curves/combined_df.Rdata\n\n",
    "Options:\n",
    "  --deep                Run all.equal() snippets for non-identical objects.\n",
    "  --full                Use all rows for distribution analysis (no sampling).\n",
    "  --sample-rows=N       Rows sampled per data frame for analysis (default: 200000).\n",
    "  --ks-sample=N         Max values per side for KS test per numeric column (default: 50000).\n",
    "  --top=N               Number of top-shift columns to print (default: 15).\n",
    "  --top-levels=N        Changed categorical levels to print per column (default: 5).\n",
    "  --object=NAME         Restrict comparison to a single object name.\n",
    "  --help                Show this message.\n",
    sep = ""
  )
}

parse_int_opt <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 0L) {
    stop("Invalid value for ", name, ": ", value)
  }
  out
}

opts <- list(
  deep = FALSE,
  sample_rows = 200000L,
  ks_sample = 50000L,
  top = 15L,
  top_levels = 5L,
  object = NULL
)

positional <- character(0)
for (arg in args) {
  if (arg == "--help") {
    print_usage()
    quit(status = 0)
  } else if (arg == "--deep") {
    opts$deep <- TRUE
  } else if (arg == "--full") {
    opts$sample_rows <- 0L
  } else if (startsWith(arg, "--sample-rows=")) {
    opts$sample_rows <- parse_int_opt(sub("^--sample-rows=", "", arg), "--sample-rows")
  } else if (startsWith(arg, "--ks-sample=")) {
    opts$ks_sample <- parse_int_opt(sub("^--ks-sample=", "", arg), "--ks-sample")
  } else if (startsWith(arg, "--top=")) {
    opts$top <- parse_int_opt(sub("^--top=", "", arg), "--top")
  } else if (startsWith(arg, "--top-levels=")) {
    opts$top_levels <- parse_int_opt(sub("^--top-levels=", "", arg), "--top-levels")
  } else if (startsWith(arg, "--object=")) {
    opts$object <- sub("^--object=", "", arg)
    if (!nzchar(opts$object)) stop("--object must not be empty")
  } else if (startsWith(arg, "--")) {
    stop("Unknown option: ", arg)
  } else {
    positional <- c(positional, arg)
  }
}

default_a <- "data/calibration/carbon_capture_curves/combined_df.Rdata"
default_b <- "data/backup/calibration/carbon_capture_curves/combined_df.Rdata"

file_a <- if (length(positional) >= 1) positional[[1]] else default_a
file_b <- if (length(positional) >= 2) positional[[2]] else default_b

if (!file.exists(file_a)) stop("File not found: ", file_a)
if (!file.exists(file_b)) stop("File not found: ", file_b)

fmt_bytes <- function(x) {
  if (is.na(x)) return("NA")
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1L
  y <- as.numeric(x)
  while (y >= 1024 && i < length(units)) {
    y <- y / 1024
    i <- i + 1L
  }
  sprintf("%.2f %s", y, units[[i]])
}

obj_dim <- function(x) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if (is.null(d)) "NULL" else paste(d, collapse = "x")
}

obj_summary <- function(x) {
  paste0(
    "class=", paste(class(x), collapse = "/"),
    ", typeof=", typeof(x),
    ", length=", length(x),
    ", dim=", obj_dim(x)
  )
}

sample_idx <- function(n, sample_rows) {
  if (sample_rows <= 0L || n <= sample_rows) return(seq_len(n))
  sort(sample.int(n, sample_rows, replace = FALSE))
}

is_numeric_like <- function(x) {
  is.numeric(x) || inherits(x, "Date") || inherits(x, "POSIXt")
}

as_numeric_clean <- function(x) {
  as.numeric(x)
}

numeric_stats <- function(x) {
  x <- as_numeric_clean(x)
  n_all <- length(x)
  n_na <- sum(is.na(x))
  x <- x[!is.na(x)]

  out <- list(
    n_all = n_all,
    n_non_na = length(x),
    na_rate = if (n_all == 0) NA_real_ else n_na / n_all,
    mean = NA_real_,
    sd = NA_real_,
    min = NA_real_,
    q01 = NA_real_,
    q05 = NA_real_,
    q25 = NA_real_,
    q50 = NA_real_,
    q75 = NA_real_,
    q95 = NA_real_,
    q99 = NA_real_,
    max = NA_real_
  )

  if (length(x) == 0L) return(out)

  q <- as.numeric(stats::quantile(x, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99), na.rm = TRUE, names = FALSE, type = 7))
  out$mean <- mean(x)
  out$sd <- stats::sd(x)
  out$min <- min(x)
  out$q01 <- q[[1]]
  out$q05 <- q[[2]]
  out$q25 <- q[[3]]
  out$q50 <- q[[4]]
  out$q75 <- q[[5]]
  out$q95 <- q[[6]]
  out$q99 <- q[[7]]
  out$max <- max(x)
  out
}

safe_ks <- function(x, y, ks_sample) {
  x <- as_numeric_clean(x)
  y <- as_numeric_clean(y)
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2L || length(y) < 2L) {
    return(c(stat = NA_real_, p = NA_real_))
  }
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(c(stat = NA_real_, p = NA_real_))
  }
  if (ks_sample > 0L && length(x) > ks_sample) x <- sample(x, ks_sample)
  if (ks_sample > 0L && length(y) > ks_sample) y <- sample(y, ks_sample)

  ks <- suppressWarnings(tryCatch(stats::ks.test(x, y), error = function(e) NULL))
  if (is.null(ks)) return(c(stat = NA_real_, p = NA_real_))
  c(stat = unname(ks$statistic), p = ks$p.value)
}

format_num <- function(x) {
  if (is.na(x)) return("NA")
  formatC(x, digits = 6, format = "fg", flag = "#")
}

categorical_profile <- function(x) {
  y <- if (is.factor(x)) as.character(x) else as.character(x)
  y[is.na(y)] <- "<NA>"
  tab <- table(y, useNA = "no")
  prop <- tab / sum(tab)
  list(prop = prop, distinct = length(tab), na_rate = mean(is.na(x)))
}

compare_categorical <- function(x, y, top_levels) {
  px <- categorical_profile(x)
  py <- categorical_profile(y)

  lvls <- union(names(px$prop), names(py$prop))
  ax <- setNames(rep(0, length(lvls)), lvls)
  by <- setNames(rep(0, length(lvls)), lvls)
  ax[names(px$prop)] <- as.numeric(px$prop)
  by[names(py$prop)] <- as.numeric(py$prop)

  delta <- ax - by
  tvd <- 0.5 * sum(abs(delta))
  ord <- order(abs(delta), decreasing = TRUE)
  keep <- seq_len(min(length(ord), max(1L, top_levels)))

  list(
    tvd = tvd,
    distinct_a = px$distinct,
    distinct_b = py$distinct,
    na_rate_a = px$na_rate,
    na_rate_b = py$na_rate,
    max_level_shift = if (length(delta) == 0L) 0 else max(abs(delta)),
    top_levels = data.frame(
      level = names(delta)[ord][keep],
      delta_share = as.numeric(delta[ord][keep]),
      stringsAsFactors = FALSE
    )
  )
}

compare_data_frames <- function(a, b, object_name, opts) {
  cat("  data.frame_structure:\n")
  cat("    nrow_a=", nrow(a), ", ncol_a=", ncol(a), "\n", sep = "")
  cat("    nrow_b=", nrow(b), ", ncol_b=", ncol(b), "\n", sep = "")

  cols_a <- names(a)
  cols_b <- names(b)
  only_a <- setdiff(cols_a, cols_b)
  only_b <- setdiff(cols_b, cols_a)
  common <- intersect(cols_a, cols_b)

  cat("    cols_only_in_a=", if (length(only_a) == 0L) "<none>" else paste(only_a, collapse = ", "), "\n", sep = "")
  cat("    cols_only_in_b=", if (length(only_b) == 0L) "<none>" else paste(only_b, collapse = ", "), "\n", sep = "")
  cat("    shared_cols=", length(common), "\n", sep = "")

  if (length(common) == 0L) {
    cat("  no shared columns; skipping distribution comparison\n")
    return(invisible(NULL))
  }

  idx_a <- sample_idx(nrow(a), opts$sample_rows)
  idx_b <- sample_idx(nrow(b), opts$sample_rows)
  a_s <- a[idx_a, common, drop = FALSE]
  b_s <- b[idx_b, common, drop = FALSE]

  cat("  sampled_rows_a=", nrow(a_s), " / ", nrow(a), "\n", sep = "")
  cat("  sampled_rows_b=", nrow(b_s), " / ", nrow(b), "\n", sep = "")

  class_changed <- vapply(
    common,
    function(cl) paste(class(a[[cl]]), collapse = "/") != paste(class(b[[cl]]), collapse = "/"),
    logical(1)
  )
  changed_cols <- common[class_changed]
  cat("  shared_cols_with_class_change=", length(changed_cols), "\n", sep = "")
  if (length(changed_cols) > 0L) {
    show <- head(changed_cols, opts$top)
    for (cl in show) {
      cat("    class_change ", cl, ": ", paste(class(a[[cl]]), collapse = "/"), " -> ", paste(class(b[[cl]]), collapse = "/"), "\n", sep = "")
    }
  }

  numeric_cols <- common[vapply(common, function(cl) is_numeric_like(a_s[[cl]]) && is_numeric_like(b_s[[cl]]), logical(1))]
  categorical_cols <- setdiff(common, numeric_cols)

  cat("  numeric_shared_cols=", length(numeric_cols), "\n", sep = "")
  cat("  categorical_shared_cols=", length(categorical_cols), "\n", sep = "")

  if (length(numeric_cols) > 0L) {
    num_rows <- vector("list", length(numeric_cols))

    for (i in seq_along(numeric_cols)) {
      cl <- numeric_cols[[i]]
      sa <- numeric_stats(a_s[[cl]])
      sb <- numeric_stats(b_s[[cl]])
      ks <- safe_ks(a_s[[cl]], b_s[[cl]], opts$ks_sample)

      pooled_sd <- sqrt((sa$sd ^ 2 + sb$sd ^ 2) / 2)
      cohen_d <- if (is.na(pooled_sd) || pooled_sd == 0) NA_real_ else abs(sa$mean - sb$mean) / pooled_sd
      iqr_b <- sb$q75 - sb$q25
      median_shift_iqr <- if (is.na(iqr_b) || abs(iqr_b) < 1e-12) NA_real_ else abs(sa$q50 - sb$q50) / abs(iqr_b)
      na_diff <- abs(sa$na_rate - sb$na_rate)
      score <- suppressWarnings(max(c(cohen_d, median_shift_iqr, na_diff * 5), na.rm = TRUE))
      if (!is.finite(score)) score <- 0

      num_rows[[i]] <- data.frame(
        column = cl,
        score = score,
        mean_a = sa$mean,
        mean_b = sb$mean,
        median_a = sa$q50,
        median_b = sb$q50,
        iqr_a = sa$q75 - sa$q25,
        iqr_b = sb$q75 - sb$q25,
        p95_a = sa$q95,
        p95_b = sb$q95,
        na_rate_a = sa$na_rate,
        na_rate_b = sb$na_rate,
        cohen_d = cohen_d,
        median_shift_iqr = median_shift_iqr,
        ks_stat = ks[["stat"]],
        ks_p = ks[["p"]],
        stringsAsFactors = FALSE
      )
    }

    num_df <- do.call(rbind, num_rows)
    num_df <- num_df[order(num_df$score, decreasing = TRUE), , drop = FALSE]

    cat("\n  numeric_distribution_shift_top=", min(nrow(num_df), opts$top), "\n", sep = "")
    show_n <- seq_len(min(nrow(num_df), opts$top))
    for (i in show_n) {
      r <- num_df[i, ]
      cat(
        "    ", r$column,
        " | score=", format_num(r$score),
        " | mean ", format_num(r$mean_a), " -> ", format_num(r$mean_b),
        " | median ", format_num(r$median_a), " -> ", format_num(r$median_b),
        " | IQR ", format_num(r$iqr_a), " -> ", format_num(r$iqr_b),
        " | p95 ", format_num(r$p95_a), " -> ", format_num(r$p95_b),
        " | NA ", format_num(r$na_rate_a), " -> ", format_num(r$na_rate_b),
        " | d=", format_num(r$cohen_d),
        " | med/IQR=", format_num(r$median_shift_iqr),
        " | KS(stat,p)=(", format_num(r$ks_stat), ", ", format_num(r$ks_p), ")",
        "\n",
        sep = ""
      )
    }
  }

  if (length(categorical_cols) > 0L) {
    cat_rows <- vector("list", length(categorical_cols))
    cat_detail <- vector("list", length(categorical_cols))

    for (i in seq_along(categorical_cols)) {
      cl <- categorical_cols[[i]]
      cmp <- compare_categorical(a_s[[cl]], b_s[[cl]], opts$top_levels)
      cat_rows[[i]] <- data.frame(
        column = cl,
        tvd = cmp$tvd,
        distinct_a = cmp$distinct_a,
        distinct_b = cmp$distinct_b,
        na_rate_a = cmp$na_rate_a,
        na_rate_b = cmp$na_rate_b,
        max_level_shift = cmp$max_level_shift,
        stringsAsFactors = FALSE
      )
      cat_detail[[i]] <- cmp$top_levels
    }

    cat_df <- do.call(rbind, cat_rows)
    cat_df <- cat_df[order(cat_df$tvd, decreasing = TRUE), , drop = FALSE]

    cat("\n  categorical_distribution_shift_top=", min(nrow(cat_df), opts$top), "\n", sep = "")
    show_n <- seq_len(min(nrow(cat_df), opts$top))
    for (i in show_n) {
      r <- cat_df[i, ]
      cat(
        "    ", r$column,
        " | TVD=", format_num(r$tvd),
        " | distinct ", r$distinct_a, " -> ", r$distinct_b,
        " | NA ", format_num(r$na_rate_a), " -> ", format_num(r$na_rate_b),
        " | max_level_shift=", format_num(r$max_level_shift),
        "\n",
        sep = ""
      )

      idx <- match(r$column, categorical_cols)
      lev <- cat_detail[[idx]]
      if (!is.null(lev) && nrow(lev) > 0L) {
        for (j in seq_len(nrow(lev))) {
          cat("      level ", lev$level[j], " delta_share=", format_num(lev$delta_share[j]), "\n", sep = "")
        }
      }
    }
  }

  if (opts$deep) {
    ae <- tryCatch(all.equal(a, b, check.attributes = TRUE), error = function(e) e$message)
    if (isTRUE(ae)) {
      cat("\n  all.equal: TRUE\n")
    } else if (is.character(ae)) {
      n <- min(length(ae), 10L)
      cat("\n  all.equal (first ", n, "): ", paste(ae[seq_len(n)], collapse = " | "), "\n", sep = "")
    } else {
      cat("\n  all.equal: ", as.character(ae), "\n", sep = "")
    }
  }

  invisible(NULL)
}

cat("== Comparison Setup ==\n")
cat("file_a=", file_a, "\n", sep = "")
cat("file_b=", file_b, "\n", sep = "")
cat("options: deep=", opts$deep,
    ", sample_rows=", opts$sample_rows,
    ", ks_sample=", opts$ks_sample,
    ", top=", opts$top,
    ", top_levels=", opts$top_levels,
    if (!is.null(opts$object)) paste0(", object=", opts$object) else "",
    "\n\n", sep = "")

cat("== File-level comparison ==\n")
fi_a <- file.info(file_a)
fi_b <- file.info(file_b)
md5 <- tools::md5sum(c(file_a, file_b))

cat("A: ", file_a, "\n", sep = "")
cat("   size=", fi_a$size, " (", fmt_bytes(fi_a$size), ")\n", sep = "")
cat("   mtime=", as.character(fi_a$mtime), "\n", sep = "")
cat("   md5=", unname(md5[[file_a]]), "\n", sep = "")

cat("B: ", file_b, "\n", sep = "")
cat("   size=", fi_b$size, " (", fmt_bytes(fi_b$size), ")\n", sep = "")
cat("   mtime=", as.character(fi_b$mtime), "\n", sep = "")
cat("   md5=", unname(md5[[file_b]]), "\n", sep = "")

same_file <- identical(unname(md5[[file_a]]), unname(md5[[file_b]]))
cat("same_bytes=", same_file, "\n\n", sep = "")

cat("== Load .Rdata files ==\n")
env_a <- new.env(parent = emptyenv())
env_b <- new.env(parent = emptyenv())

time_a <- system.time(load(file_a, envir = env_a))
time_b <- system.time(load(file_b, envir = env_b))

cat("load_time_a_sec=", unname(time_a[["elapsed"]]), "\n", sep = "")
cat("load_time_b_sec=", unname(time_b[["elapsed"]]), "\n\n", sep = "")

names_a <- sort(ls(env_a, all.names = TRUE))
names_b <- sort(ls(env_b, all.names = TRUE))

cat("== Object inventory ==\n")
cat("objects_in_a=", length(names_a), "\n", sep = "")
cat("objects_in_b=", length(names_b), "\n", sep = "")

only_a <- setdiff(names_a, names_b)
only_b <- setdiff(names_b, names_a)
common <- intersect(names_a, names_b)

cat("only_in_a=", if (length(only_a) == 0L) "<none>" else paste(only_a, collapse = ", "), "\n", sep = "")
cat("only_in_b=", if (length(only_b) == 0L) "<none>" else paste(only_b, collapse = ", "), "\n", sep = "")

if (!is.null(opts$object)) {
  common <- intersect(common, opts$object)
  if (length(common) == 0L) {
    stop("Requested object not found in both files: ", opts$object)
  }
}

cat("common_objects_to_compare=", length(common), "\n\n", sep = "")

if (length(common) == 0L) {
  cat("No common objects to compare.\n")
  quit(status = 0)
}

different_count <- 0L

for (nm in common) {
  cat("== Object: ", nm, " ==\n", sep = "")
  a <- get(nm, envir = env_a, inherits = FALSE)
  b <- get(nm, envir = env_b, inherits = FALSE)

  same <- identical(a, b)
  cat("identical=", same, "\n", sep = "")
  cat("A: ", obj_summary(a), "\n", sep = "")
  cat("B: ", obj_summary(b), "\n", sep = "")

  if (!same) {
    different_count <- different_count + 1L

    if (inherits(a, "data.frame") && inherits(b, "data.frame")) {
      compare_data_frames(a, b, nm, opts)
    } else {
      if (is.atomic(a) && is.atomic(b) && is_numeric_like(a) && is_numeric_like(b)) {
        sa <- numeric_stats(a)
        sb <- numeric_stats(b)
        ks <- safe_ks(a, b, opts$ks_sample)
        cat("  numeric_vector_summary:\n")
        cat("    mean ", format_num(sa$mean), " -> ", format_num(sb$mean), "\n", sep = "")
        cat("    median ", format_num(sa$q50), " -> ", format_num(sb$q50), "\n", sep = "")
        cat("    p95 ", format_num(sa$q95), " -> ", format_num(sb$q95), "\n", sep = "")
        cat("    NA ", format_num(sa$na_rate), " -> ", format_num(sb$na_rate), "\n", sep = "")
        cat("    KS(stat,p)=(", format_num(ks[["stat"]]), ", ", format_num(ks[["p"]]), ")\n", sep = "")
      } else if (is.atomic(a) && is.atomic(b)) {
        cmp <- compare_categorical(a, b, opts$top_levels)
        cat("  atomic_categorical_summary:\n")
        cat("    TVD=", format_num(cmp$tvd), "\n", sep = "")
        cat("    distinct ", cmp$distinct_a, " -> ", cmp$distinct_b, "\n", sep = "")
        cat("    NA ", format_num(cmp$na_rate_a), " -> ", format_num(cmp$na_rate_b), "\n", sep = "")
        if (!is.null(cmp$top_levels) && nrow(cmp$top_levels) > 0L) {
          for (i in seq_len(nrow(cmp$top_levels))) {
            cat("      level ", cmp$top_levels$level[i], " delta_share=", format_num(cmp$top_levels$delta_share[i]), "\n", sep = "")
          }
        }
      } else {
        cat("  non-data-frame object differs; use --deep for all.equal details.\n")
      }

      if (opts$deep) {
        ae <- tryCatch(all.equal(a, b, check.attributes = TRUE), error = function(e) e$message)
        if (isTRUE(ae)) {
          cat("  all.equal: TRUE\n")
        } else if (is.character(ae)) {
          n <- min(length(ae), 10L)
          cat("  all.equal (first ", n, "): ", paste(ae[seq_len(n)], collapse = " | "), "\n", sep = "")
        } else {
          cat("  all.equal: ", as.character(ae), "\n", sep = "")
        }
      }
    }
  }

  cat("\n")
}

cat("== Summary ==\n")
cat("different_common_objects=", different_count, " / ", length(common), "\n", sep = "")
if (length(only_a) > 0L || length(only_b) > 0L || different_count > 0L) {
  cat("RESULT: files differ\n")
} else {
  cat("RESULT: files are equivalent at object level\n")
}
