#!/usr/bin/env Rscript

# Stress-responsive dark-gene plots
# ---------------------------------
# Produces a second, publication-facing plot set focused on RNA-seq responsive
# dark genes after annotation, BUSCO, repeat/TE, DMR and DE integration.
#
# The script uses ggplot2 and applies bbplot::bbc_style() when available.
# If bbplot is not installed, it falls back to a clean minimal theme.
#
# Layout note:
# bbplot/BBC-style themes use large typography by default. The final theme layer
# below deliberately scales text, labels and margins down so long biological
# category labels fit cleanly in the exported plotting space.

options(stringsAsFactors = FALSE)
options(warn = 1)

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

if (has_pkg("readr")) {
    read_tsv_flex <- function(path) as.data.frame(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
} else {
    read_tsv_flex <- function(path) read.delim(path, sep = "\t", header = TRUE, quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
}

if (has_pkg("dplyr")) {
    library(dplyr)
}
if (has_pkg("tidyr")) {
    library(tidyr)
}
if (has_pkg("stringr")) {
    library(stringr)
}

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/uoa/scratch/users/r02hw22/project_dark_genes")
DE_CONTEXT_TSV <- Sys.getenv("DE_CONTEXT_TSV", file.path(PROJECT_DIR, "11_expression_context/equina_dark_candidates.de_context.tsv"))
DE_SIG_LONG_TSV <- Sys.getenv("DE_SIG_LONG_TSV", file.path(PROJECT_DIR, "11_expression_context/equina_dark_candidates.de_significant_long.tsv"))
OUTDIR <- Sys.getenv("STRESS_FIGURE_DIR", file.path(PROJECT_DIR, "12_final_candidates/stress_responsive_plots"))
TABLE_DIR <- file.path(OUTDIR, "plot_tables")

PADJ_THRESHOLD <- as.numeric(Sys.getenv("PADJ_THRESHOLD", "0.05"))
LFC_THRESHOLD <- as.numeric(Sys.getenv("LFC_THRESHOLD", "1"))

# Main contrasts to highlight in diesel/stress-responsive interpretation.
FOCAL_CONTRASTS <- strsplit(Sys.getenv(
    "FOCAL_CONTRASTS",
    "diesel_added_wald,combined_wald,full_model_LRT,interactive_only_wald,interactive_only_LRT"
), ",")[[1]]
FOCAL_CONTRASTS <- trimws(FOCAL_CONTRASTS)

# Contrasts shown in ordering for figures when present.
CONTRAST_ORDER <- c(
    "diesel_added_wald",
    "combined_wald",
    "full_model_LRT",
    "interactive_only_wald",
    "interactive_only_LRT",
    "diesel_only_wald",
    "diesel_only_LRT",
    "salinity_added_wald",
    "salinity_only_wald",
    "salinity_only_LRT"
)

if (!file.exists(DE_CONTEXT_TSV)) stop("DE context table not found: ", DE_CONTEXT_TSV)
if (!file.exists(DE_SIG_LONG_TSV)) stop("Significant DE long table not found: ", DE_SIG_LONG_TSV)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

message("DE context table: ", DE_CONTEXT_TSV)
message("DE significant long table: ", DE_SIG_LONG_TSV)
message("Output directory: ", OUTDIR)

clean_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == "." | x == "NA" | x == "none"] <- "NA"
    x
}

as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))

pretty_contrast <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("wald", "Wald", x, ignore.case = TRUE)
    x <- gsub("lrt", "LRT", x, ignore.case = TRUE)
    x <- gsub("diesel", "Diesel", x, ignore.case = TRUE)
    x <- gsub("salinity", "Salinity", x, ignore.case = TRUE)
    x <- gsub("combined", "Combined", x, ignore.case = TRUE)
    x <- gsub("interactive", "Interaction", x, ignore.case = TRUE)
    x
}

wrap_label <- function(x, width = 28) {
    x <- as.character(x)
    if (has_pkg("stringr")) return(stringr::str_wrap(x, width = width))
    vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

short_signature <- function(x) {
    x <- gsub("diesel_added_wald", "Diesel added", x)
    x <- gsub("combined_wald", "Combined", x)
    x <- gsub("full_model_LRT", "Full model", x)
    x <- gsub("interactive_only_wald", "Interaction Wald", x)
    x <- gsub("interactive_only_LRT", "Interaction LRT", x)
    x <- gsub("_", " ", x)
    x
}

safe_col <- function(df, candidates, default = NA_character_) {
    nms <- names(df)
    lower <- tolower(nms)
    for (cand in candidates) {
        hit <- which(lower == tolower(cand))
        if (length(hit)) return(nms[hit[1]])
    }
    default
}

style_plot <- function(p, base_size = 10) {
    # bbplot package from BBC uses bbc_style(); some users refer to this as bbplot styling.
    if (has_pkg("bbplot") && "bbc_style" %in% getNamespaceExports("bbplot")) {
        p <- p + bbplot::bbc_style()
    } else if (has_pkg("bbplot") && "bbplot" %in% getNamespaceExports("bbplot")) {
        # Defensive fallback if a local bbplot() wrapper exists in a site package.
        p <- bbplot::bbplot(p)
    } else {
        p <- p + theme_minimal(base_size = base_size)
    }

    # Applied after bbplot so oversized defaults are scaled down consistently.
    p + theme(
        plot.title = element_text(face = "bold", size = base_size + 5, lineheight = 0.95, margin = margin(b = 5)),
        plot.subtitle = element_text(size = base_size + 1, colour = "grey30", lineheight = 0.95, margin = margin(b = 8)),
        plot.caption = element_text(size = base_size - 2, colour = "grey40", hjust = 0, lineheight = 0.95),
        axis.title = element_text(face = "bold", size = base_size + 1),
        axis.text = element_text(colour = "grey15", size = base_size, lineheight = 0.9),
        axis.text.y = element_text(size = base_size - 1, lineheight = 0.9),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = base_size - 1),
        legend.key.size = unit(0.48, "lines"),
        strip.text = element_text(face = "bold", size = base_size),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 16, 10, 10)
    )
}

save_plot <- function(plot, stub, width = 11.5, height = 7.2) {
    png_path <- file.path(OUTDIR, paste0(stub, ".png"))
    pdf_path <- file.path(OUTDIR, paste0(stub, ".pdf"))
    ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
    ggsave(pdf_path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
    data.frame(stub = stub, png = png_path, pdf = pdf_path, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

contrast_factor <- function(x) {
    present_order <- c(CONTRAST_ORDER[CONTRAST_ORDER %in% unique(x)], sort(setdiff(unique(x), CONTRAST_ORDER)))
    factor(x, levels = present_order)
}

context <- read_tsv_flex(DE_CONTEXT_TSV)
sig <- read_tsv_flex(DE_SIG_LONG_TSV)

# Normalise key columns.
if (!"contrast" %in% names(sig)) stop("Significant long table must contain a contrast column")
if (!"direction" %in% names(sig)) stop("Significant long table must contain a direction column")
if (!"logfc" %in% names(sig)) stop("Significant long table must contain a logfc column")
if (!"padj" %in% names(sig)) stop("Significant long table must contain a padj column")

sig$contrast <- clean_chr(sig$contrast)
sig$direction <- clean_chr(sig$direction)
sig$logfc <- as_num(sig$logfc)
sig$padj <- as_num(sig$padj)
sig$contrast_pretty <- pretty_contrast(sig$contrast)
sig$contrast_pretty_wrapped <- wrap_label(sig$contrast_pretty, 22)
sig$contrast_factor <- contrast_factor(sig$contrast)

candidate_col <- safe_col(context, c("protein_id", "candidate_id", "transcript_id", "gene_id"), names(context)[1])
priority_col <- safe_col(context, c("priority_tier"), NA_character_)
dup_col <- safe_col(context, c("duplication_interpretation", "duplication_status"), NA_character_)
repeat_col <- safe_col(context, c("repeat_overlap_status"), NA_character_)
dmr_col <- safe_col(context, c("dmr_overlap_status"), NA_character_)
de_status_col <- safe_col(context, c("de_match_status"), NA_character_)

context$candidate_id_for_plot <- clean_chr(context[[candidate_col]])
if (!is.na(priority_col)) context$priority_tier_for_plot <- clean_chr(context[[priority_col]]) else context$priority_tier_for_plot <- "unknown_priority"
if (!is.na(dup_col)) context$duplication_for_plot <- clean_chr(context[[dup_col]]) else context$duplication_for_plot <- "unknown_duplication"
if (!is.na(repeat_col)) context$repeat_for_plot <- clean_chr(context[[repeat_col]]) else context$repeat_for_plot <- "unknown_repeat_status"
if (!is.na(dmr_col)) context$dmr_for_plot <- clean_chr(context[[dmr_col]]) else context$dmr_for_plot <- "unknown_dmr_status"
if (!is.na(de_status_col)) context$de_status_for_plot <- clean_chr(context[[de_status_col]]) else context$de_status_for_plot <- "unknown_de_status"

context$priority_tier_wrapped <- wrap_label(gsub("_", " ", context$priority_tier_for_plot), 24)
context$repeat_status_wrapped <- wrap_label(gsub("_", " ", context$repeat_for_plot), 24)

sig_candidates <- unique(clean_chr(sig$candidate_id))
context$stress_responsive_any <- context$candidate_id_for_plot %in% sig_candidates

focal_sig <- sig[sig$contrast %in% FOCAL_CONTRASTS, , drop = FALSE]
focal_candidates <- unique(clean_chr(focal_sig$candidate_id))
context$stress_responsive_focal <- context$candidate_id_for_plot %in% focal_candidates

# Build compact plot tables.
contrast_counts <- aggregate(candidate_id ~ contrast + contrast_pretty, unique(sig[, c("candidate_id", "contrast", "contrast_pretty")]), length)
names(contrast_counts)[names(contrast_counts) == "candidate_id"] <- "n_dark_candidates"
contrast_counts$contrast_factor <- contrast_factor(contrast_counts$contrast)
contrast_counts$contrast_pretty_wrapped <- wrap_label(contrast_counts$contrast_pretty, 22)
contrast_counts <- contrast_counts[order(contrast_counts$contrast_factor), ]
write_tsv(contrast_counts[, c("contrast", "contrast_pretty", "n_dark_candidates")], file.path(TABLE_DIR, "stress_figure_16_de_candidates_by_contrast.tsv"))

direction_counts <- aggregate(candidate_id ~ contrast + contrast_pretty + direction, unique(sig[, c("candidate_id", "contrast", "contrast_pretty", "direction")]), length)
names(direction_counts)[names(direction_counts) == "candidate_id"] <- "n_dark_candidates"
direction_counts$contrast_factor <- contrast_factor(direction_counts$contrast)
direction_counts$contrast_pretty_wrapped <- wrap_label(direction_counts$contrast_pretty, 22)
write_tsv(direction_counts[, c("contrast", "contrast_pretty", "direction", "n_dark_candidates")], file.path(TABLE_DIR, "stress_figure_17_direction_by_contrast.tsv"))

priority_stress <- as.data.frame(table(context$priority_tier_wrapped, context$stress_responsive_focal), stringsAsFactors = FALSE)
names(priority_stress) <- c("priority_tier", "focal_stress_responsive", "n_dark_candidates")
priority_stress$focal_stress_responsive <- ifelse(priority_stress$focal_stress_responsive == "TRUE", "Focal stress-responsive", "Not focal stress-responsive")
write_tsv(priority_stress, file.path(TABLE_DIR, "stress_figure_18_priority_by_focal_stress.tsv"))

repeat_stress <- as.data.frame(table(context$repeat_status_wrapped, context$stress_responsive_focal), stringsAsFactors = FALSE)
names(repeat_stress) <- c("repeat_status", "focal_stress_responsive", "n_dark_candidates")
repeat_stress$focal_stress_responsive <- ifelse(repeat_stress$focal_stress_responsive == "TRUE", "Focal stress-responsive", "Not focal stress-responsive")
write_tsv(repeat_stress, file.path(TABLE_DIR, "stress_figure_19_repeat_by_focal_stress.tsv"))

# Multi-contrast membership table for focal contrasts.
membership <- unique(sig[, c("candidate_id", "contrast")])
membership <- membership[membership$contrast %in% FOCAL_CONTRASTS, , drop = FALSE]
if (nrow(membership) > 0) {
    membership$present <- 1L
    if (has_pkg("tidyr")) {
        membership_wide <- tidyr::pivot_wider(membership, names_from = contrast, values_from = present, values_fill = 0)
        membership_wide <- as.data.frame(membership_wide)
    } else {
        candidate_ids <- unique(membership$candidate_id)
        membership_wide <- data.frame(candidate_id = candidate_ids, stringsAsFactors = FALSE)
        for (ct in unique(membership$contrast)) {
            membership_wide[[ct]] <- as.integer(candidate_ids %in% membership$candidate_id[membership$contrast == ct])
        }
    }
    contrast_cols <- setdiff(names(membership_wide), "candidate_id")
    membership_wide$signature <- apply(membership_wide[, contrast_cols, drop = FALSE], 1, function(z) paste(contrast_cols[as.integer(z) == 1], collapse = " + "))
    membership_wide$signature_pretty <- wrap_label(short_signature(membership_wide$signature), 34)
    signature_counts <- as.data.frame(table(membership_wide$signature, membership_wide$signature_pretty), stringsAsFactors = FALSE)
    names(signature_counts) <- c("stress_signature", "stress_signature_pretty", "n_dark_candidates")
    signature_counts <- signature_counts[order(signature_counts$n_dark_candidates, decreasing = TRUE), ]
} else {
    signature_counts <- data.frame(stress_signature = character(), stress_signature_pretty = character(), n_dark_candidates = integer())
}
write_tsv(signature_counts, file.path(TABLE_DIR, "stress_figure_20_focal_stress_signature_counts.tsv"))

# Diesel-specific volcano/significance table from significant long table only.
diesel_sig <- sig[sig$contrast == "diesel_added_wald", , drop = FALSE]
diesel_sig$minus_log10_padj <- -log10(pmax(diesel_sig$padj, .Machine$double.xmin))
diesel_sig$label <- diesel_sig$candidate_id
write_tsv(diesel_sig[, c("candidate_id", "matched_de_id", "direction", "logfc", "padj", "minus_log10_padj")], file.path(TABLE_DIR, "stress_figure_21_diesel_added_significant_candidates.tsv"))

# Plot 16: significant dark candidates by contrast.
p16 <- ggplot(contrast_counts, aes(x = reorder(contrast_pretty_wrapped, n_dark_candidates), y = n_dark_candidates)) +
    geom_col(width = 0.60) +
    coord_flip(clip = "off") +
    geom_text(aes(label = n_dark_candidates), hjust = -0.10, size = 2.8) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        title = "Dark genes respond across the multi-stressor RNA-seq design",
        subtitle = paste0("Significant when padj ≤ ", PADJ_THRESHOLD, " and |logFC| ≥ ", LFC_THRESHOLD),
        x = NULL,
        y = "Significant dark candidates",
        caption = "Source: integrated dark-candidate DE context table"
    )
p16 <- style_plot(p16, base_size = 9)

# Plot 17: direction of response by contrast.
p17 <- ggplot(direction_counts, aes(x = contrast_pretty_wrapped, y = n_dark_candidates, fill = direction)) +
    geom_col(width = 0.60, position = "stack") +
    coord_flip(clip = "off") +
    labs(
        title = "Direction of differential expression among dark genes",
        subtitle = "Up- and down-regulated significant dark candidates by contrast",
        x = NULL,
        y = "Significant dark candidates",
        caption = "Direction is based on the sign of logFC"
    )
p17 <- style_plot(p17, base_size = 9)

# Plot 18: focal stress responsiveness across priority tiers.
p18 <- ggplot(priority_stress, aes(x = priority_tier, y = n_dark_candidates, fill = focal_stress_responsive)) +
    geom_col(width = 0.58) +
    coord_flip(clip = "off") +
    labs(
        title = "Stress-responsive dark genes by priority tier",
        subtitle = "Focal stress contrasts include diesel-added, combined, full-model and interaction contrasts",
        x = NULL,
        y = "Dark candidates",
        caption = "Priority tiers are BUSCO-backed and candidate-level"
    )
p18 <- style_plot(p18, base_size = 9)

# Plot 19: repeat/TE context of focal stress-responsive dark genes.
p19 <- ggplot(repeat_stress, aes(x = repeat_status, y = n_dark_candidates, fill = focal_stress_responsive)) +
    geom_col(width = 0.58) +
    coord_flip(clip = "off") +
    labs(
        title = "Stress-responsive dark genes are shown with repeat/TE context",
        subtitle = "Repeat overlap is retained as biological context, not as an automatic exclusion",
        x = NULL,
        y = "Dark candidates",
        caption = "Repeat/TE overlap used filtered repeat-like GFF3 feature types"
    )
p19 <- style_plot(p19, base_size = 9)

# Plot 20: focal stress signature combinations.
if (nrow(signature_counts) > 0) {
    top_signatures <- head(signature_counts, 12)
    p20 <- ggplot(top_signatures, aes(x = reorder(stress_signature_pretty, n_dark_candidates), y = n_dark_candidates)) +
        geom_col(width = 0.58) +
        coord_flip(clip = "off") +
        geom_text(aes(label = n_dark_candidates), hjust = -0.10, size = 2.6) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
        labs(
            title = "Focal stress-response signatures among dark genes",
            subtitle = "Top combinations of diesel, combined, full-model and interaction evidence",
            x = NULL,
            y = "Dark candidates",
            caption = "Only candidates significant in at least one focal stress contrast are shown"
        )
    p20 <- style_plot(p20, base_size = 8)
} else {
    p20 <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No focal stress-responsive dark candidates detected", size = 4) +
        theme_void() +
        labs(title = "Focal stress-response signatures among dark genes")
}

# Plot 21: diesel-added significant candidates, logFC by significance.
if (nrow(diesel_sig) > 0) {
    diesel_sig$direction <- factor(diesel_sig$direction, levels = c("down", "up", "zero"))
    p21 <- ggplot(diesel_sig, aes(x = logfc, y = minus_log10_padj, fill = direction)) +
        geom_point(shape = 21, size = 2.5, alpha = 0.82, colour = "grey15", stroke = 0.22) +
        geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", linewidth = 0.30) +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", linewidth = 0.30) +
        labs(
            title = "Diesel-added response highlights a focused set of dark genes",
            subtitle = paste0(nrow(diesel_sig), " significant dark candidates in diesel_added_wald"),
            x = "logFC",
            y = expression(-log[10](adjusted~p)),
            caption = "Only significant diesel_added_wald dark candidates are plotted"
        )
    if (has_pkg("ggrepel")) {
        top_lab <- diesel_sig[order(diesel_sig$padj, -abs(diesel_sig$logfc)), , drop = FALSE]
        top_lab <- head(top_lab, 8)
        p21 <- p21 + ggrepel::geom_text_repel(
            data = top_lab,
            aes(label = candidate_id),
            size = 2.2,
            max.overlaps = Inf,
            min.segment.length = 0,
            box.padding = 0.25,
            point.padding = 0.18,
            segment.size = 0.22
        )
    }
    p21 <- style_plot(p21, base_size = 9)
} else {
    p21 <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No diesel_added_wald significant dark candidates detected", size = 4) +
        theme_void() +
        labs(title = "Diesel-added response among dark genes")
}

# Plot 22: final evidence stack summary.
evidence_summary <- data.frame(
    evidence_layer = c(
        "Dark candidates extracted",
        "Matched RNA-seq DE records",
        "DE in ≥1 contrast",
        "Diesel-added DE",
        "Repeat/TE overlap",
        "No DMR overlap"
    ),
    n_dark_candidates = c(
        nrow(context),
        sum(context$de_status_for_plot != "no_de_record_matched"),
        sum(context$de_status_for_plot == "de_significant"),
        length(unique(diesel_sig$candidate_id)),
        sum(context$repeat_for_plot == "repeat_overlap"),
        nrow(context)
    ),
    stringsAsFactors = FALSE
)
evidence_summary$evidence_layer_wrapped <- wrap_label(evidence_summary$evidence_layer, 24)
evidence_summary$evidence_layer_wrapped <- factor(evidence_summary$evidence_layer_wrapped, levels = rev(evidence_summary$evidence_layer_wrapped))
write_tsv(evidence_summary[, c("evidence_layer", "n_dark_candidates")], file.path(TABLE_DIR, "stress_figure_22_evidence_stack_summary.tsv"))

p22 <- ggplot(evidence_summary, aes(x = evidence_layer_wrapped, y = n_dark_candidates)) +
    geom_col(width = 0.60) +
    coord_flip(clip = "off") +
    geom_text(aes(label = n_dark_candidates), hjust = -0.10, size = 2.8) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        title = "Evidence stack for dark-gene prioritisation",
        subtitle = "Expression evidence provides the strongest current stress-response signal",
        x = NULL,
        y = "Dark candidates",
        caption = "DMR overlap is a valid negative result for both diesel methylation experiments"
    )
p22 <- style_plot(p22, base_size = 9)

manifest <- rbind(
    save_plot(p16, "stress_figure_16_de_candidates_by_contrast", 11.8, 7.6),
    save_plot(p17, "stress_figure_17_de_direction_by_contrast", 11.8, 7.6),
    save_plot(p18, "stress_figure_18_priority_vs_focal_stress", 11.8, 6.6),
    save_plot(p19, "stress_figure_19_repeat_context_vs_focal_stress", 11.8, 6.6),
    save_plot(p20, "stress_figure_20_focal_stress_signature_combinations", 13.2, 8.2),
    save_plot(p21, "stress_figure_21_diesel_added_dark_gene_response", 10.8, 7.4),
    save_plot(p22, "stress_figure_22_dark_gene_evidence_stack", 11.8, 6.8)
)

manifest$title <- c(
    "Significant dark candidates by DE contrast",
    "Direction of significant dark-candidate DE by contrast",
    "Focal stress responsiveness by BUSCO-backed priority tier",
    "Repeat/TE context of focal stress-responsive dark candidates",
    "Focal stress-response signature combinations",
    "Diesel-added significant dark-candidate response",
    "Integrated evidence stack for dark candidates"
)
write_tsv(manifest, file.path(OUTDIR, "stress_responsive_dark_gene_plot_manifest.tsv"))

summary_path <- file.path(OUTDIR, "stress_responsive_dark_gene_plots.summary.txt")
con <- file(summary_path, open = "wt")
writeLines("Stress-responsive dark-gene plot summary", con)
writeLines("========================================", con)
writeLines(paste0("DE context table: ", DE_CONTEXT_TSV), con)
writeLines(paste0("Significant DE table: ", DE_SIG_LONG_TSV), con)
writeLines(paste0("Output directory: ", OUTDIR), con)
writeLines(paste0("Candidates read: ", nrow(context)), con)
writeLines(paste0("Significant DE rows read: ", nrow(sig)), con)
writeLines(paste0("Unique significant dark candidates: ", length(sig_candidates)), con)
writeLines(paste0("Unique focal stress-responsive dark candidates: ", length(focal_candidates)), con)
writeLines(paste0("Diesel-added significant dark candidates: ", length(unique(diesel_sig$candidate_id))), con)
writeLines("", con)
writeLines("Layout settings:", con)
writeLines("- bbplot/BBC styling applied when available", con)
writeLines("- final theme layer scales text down for crowded biological labels", con)
writeLines("- long category labels are wrapped before plotting", con)
writeLines("- plot canvases are exported wider than the first version", con)
writeLines("", con)
writeLines("Contrast counts:", con)
for (i in seq_len(nrow(contrast_counts))) {
    writeLines(paste0("- ", contrast_counts$contrast[i], ": ", contrast_counts$n_dark_candidates[i]), con)
}
writeLines("", con)
writeLines("Plot files:", con)
for (i in seq_len(nrow(manifest))) {
    writeLines(paste0("- ", manifest$stub[i], ": ", manifest$png[i]), con)
}
close(con)

message("Done.")
message("Summary: ", summary_path)
message("Manifest: ", file.path(OUTDIR, "stress_responsive_dark_gene_plot_manifest.tsv"))
