#!/bin/bash
#SBATCH --job-name=make_figures
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/

PROJECT_DIR="${PROJECT_DIR:-/uoa/scratch/users/r02hw22/project_dark_genes}"

MASTER_TSV="${MASTER_TSV:-${PROJECT_DIR}/02_annotation/master/full/equina_representative_full.master_annotation.tsv}"
DARK_TSV="${DARK_TSV:-${PROJECT_DIR}/03_dark_candidates/equina_dark_candidates.tsv}"

FIGURE_DIR="${FIGURE_DIR:-${PROJECT_DIR}/08_figures}"
TABLE_DIR="${FIGURE_DIR}/figure_tables"
SANITISED_DIR="${TABLE_DIR}/sanitised_inputs"
RSCRIPT_PATH="${PROJECT_DIR}/scripts/make_figures_02_to_15.R"
SANITISER="${PROJECT_DIR}/scripts/sanitise_tsv_for_R.py"

mkdir -p "${PROJECT_DIR}/scripts" "${FIGURE_DIR}" "${TABLE_DIR}" "${SANITISED_DIR}"

find_latest_file() {
    local filename="$1"
    find "${PROJECT_DIR}" -type f -name "${filename}" -printf "%T@ %p\n" 2>/dev/null \
        | sort -nr \
        | head -n 1 \
        | cut -d' ' -f2-
}

GENOME_LINKED_TSV="${GENOME_LINKED_TSV:-$(find_latest_file "equina_dark_candidates.genome_linked.tsv")}"
SCAFFOLD_SUMMARY_TSV="${SCAFFOLD_SUMMARY_TSV:-$(find_latest_file "equina_dark_candidates.summary_by_scaffold.tsv")}"
BUSCO_SCAFFOLD_TSV="${BUSCO_SCAFFOLD_TSV:-$(find_latest_file "equina_scaffold_busco_duplication_summary.tsv")}"
CLUSTER_TSV="${CLUSTER_TSV:-$(find_latest_file "equina_dark_candidate_protein_clusters.tsv")}"
PRIORITISED_TSV="${PRIORITISED_TSV:-$(find_latest_file "equina_dark_candidates.prioritised.tsv")}"

BUSCO_FULL_TABLE="${BUSCO_FULL_TABLE:-}"
if [[ -z "${BUSCO_FULL_TABLE}" ]]; then
    BUSCO_FULL_TABLE="$(find "${PROJECT_DIR}" -type f \( -name "full_table*.tsv" -o -name "full_table*.txt" \) -printf "%T@ %p\n" 2>/dev/null \
        | sort -nr \
        | head -n 1 \
        | cut -d' ' -f2- || true)"
fi

if [[ -z "${BUSCO_FULL_TABLE}" ]]; then
    BUSCO_FULL_TABLE="NA"
fi

if [[ -z "${GENOME_LINKED_TSV}" ]]; then GENOME_LINKED_TSV="NA"; fi
if [[ -z "${SCAFFOLD_SUMMARY_TSV}" ]]; then SCAFFOLD_SUMMARY_TSV="NA"; fi
if [[ -z "${BUSCO_SCAFFOLD_TSV}" ]]; then BUSCO_SCAFFOLD_TSV="NA"; fi
if [[ -z "${CLUSTER_TSV}" ]]; then CLUSTER_TSV="NA"; fi
if [[ -z "${PRIORITISED_TSV}" ]]; then PRIORITISED_TSV="NA"; fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found. Load R first, then rerun." >&2
    echo "Example: module load R" >&2
    exit 1
fi

cat > "${SANITISER}" <<'PYTHON'
#!/usr/bin/env python3

import argparse
import csv
import os
import sys

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    limit = sys.maxsize
    while True:
        limit //= 10
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            continue


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    if args.input in {"", "NA"} or not os.path.exists(args.input) or os.path.getsize(args.input) == 0:
        print(f"SKIP\t{args.input}\tmissing")
        return

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    long_rows = 0
    short_rows = 0
    total_rows = 0

    with open(args.input, "r", newline="", encoding="utf-8", errors="replace") as fin, \
            open(args.output, "w", newline="", encoding="utf-8") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(
            fout,
            delimiter="\t",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )

        try:
            header = next(reader)
        except StopIteration:
            print(f"SKIP\t{args.input}\tempty")
            return

        ncol = len(header)
        writer.writerow(header)

        for row in reader:
            total_rows += 1

            if len(row) < ncol:
                row = row + [""] * (ncol - len(row))
                short_rows += 1

            elif len(row) > ncol:
                # Preserve early columns and collapse any excess fields into the final column.
                # This handles annotation text fields that contain extra literal tab characters.
                row = row[:ncol - 1] + [" | ".join(row[ncol - 1:])]
                long_rows += 1

            writer.writerow(row)

    print(
        f"OK\t{args.input}\t{args.output}\trows={total_rows}\tlong_rows_fixed={long_rows}\tshort_rows_fixed={short_rows}"
    )


if __name__ == "__main__":
    main()
PYTHON

chmod +x "${SANITISER}"

sanitise_or_na() {
    local src="$1"
    local stub="$2"
    local dest="${SANITISED_DIR}/${stub}.sanitised.tsv"

    if [[ "${src}" == "NA" || -z "${src}" || ! -s "${src}" ]]; then
        echo "NA"
    else
        python3 "${SANITISER}" --input "${src}" --output "${dest}" >&2
        echo "${dest}"
    fi
}

echo "Sanitising TSV inputs for R..."

MASTER_TSV_CLEAN="$(sanitise_or_na "${MASTER_TSV}" "master_annotation")"
DARK_TSV_CLEAN="$(sanitise_or_na "${DARK_TSV}" "dark_candidates")"
GENOME_LINKED_TSV_CLEAN="$(sanitise_or_na "${GENOME_LINKED_TSV}" "genome_linked")"
SCAFFOLD_SUMMARY_TSV_CLEAN="$(sanitise_or_na "${SCAFFOLD_SUMMARY_TSV}" "scaffold_summary")"
BUSCO_SCAFFOLD_TSV_CLEAN="$(sanitise_or_na "${BUSCO_SCAFFOLD_TSV}" "busco_scaffold_summary")"
CLUSTER_TSV_CLEAN="$(sanitise_or_na "${CLUSTER_TSV}" "dark_candidate_clusters")"
PRIORITISED_TSV_CLEAN="$(sanitise_or_na "${PRIORITISED_TSV}" "prioritised_candidates")"

export PROJECT_DIR
export FIGURE_DIR
export TABLE_DIR
export MASTER_TSV="${MASTER_TSV_CLEAN}"
export DARK_TSV="${DARK_TSV_CLEAN}"
export GENOME_LINKED_TSV="${GENOME_LINKED_TSV_CLEAN}"
export SCAFFOLD_SUMMARY_TSV="${SCAFFOLD_SUMMARY_TSV_CLEAN}"
export BUSCO_FULL_TABLE
export BUSCO_SCAFFOLD_TSV="${BUSCO_SCAFFOLD_TSV_CLEAN}"
export CLUSTER_TSV="${CLUSTER_TSV_CLEAN}"
export PRIORITISED_TSV="${PRIORITISED_TSV_CLEAN}"

cat > "${RSCRIPT_PATH}" <<'RSCRIPT'
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(warn = 1)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR")
FIGURE_DIR <- Sys.getenv("FIGURE_DIR")
TABLE_DIR <- Sys.getenv("TABLE_DIR")

MASTER_TSV <- Sys.getenv("MASTER_TSV")
DARK_TSV <- Sys.getenv("DARK_TSV")
GENOME_LINKED_TSV <- Sys.getenv("GENOME_LINKED_TSV")
SCAFFOLD_SUMMARY_TSV <- Sys.getenv("SCAFFOLD_SUMMARY_TSV")
BUSCO_FULL_TABLE <- Sys.getenv("BUSCO_FULL_TABLE")
BUSCO_SCAFFOLD_TSV <- Sys.getenv("BUSCO_SCAFFOLD_TSV")
CLUSTER_TSV <- Sys.getenv("CLUSTER_TSV")
PRIORITISED_TSV <- Sys.getenv("PRIORITISED_TSV")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
    figure = character(),
    title = character(),
    status = character(),
    png = character(),
    pdf = character(),
    notes = character(),
    stringsAsFactors = FALSE
)

add_manifest <- function(fig, title, status, png = "", pdf = "", notes = "") {
    manifest <<- rbind(
        manifest,
        data.frame(
            figure = fig,
            title = title,
            status = status,
            png = png,
            pdf = pdf,
            notes = notes,
            stringsAsFactors = FALSE
        )
    )
}

read_tsv <- function(path, label) {
    if (is.na(path) || path == "" || path == "NA" || !file.exists(path) || file.info(path)$size == 0) {
        message("Skipping read: ", label, " not found: ", path)
        return(NULL)
    }

    tryCatch(
        read.delim(
            path,
            sep = "\t",
            header = TRUE,
            quote = "\"",
            comment.char = "",
            check.names = FALSE,
            fill = TRUE,
            stringsAsFactors = FALSE
        ),
        error = function(e) {
            message("Could not read ", label, ": ", conditionMessage(e))
            NULL
        }
    )
}

detect_col <- function(df, candidates) {
    if (is.null(df)) return(NA_character_)

    nms <- names(df)
    lower <- tolower(nms)

    for (cand in candidates) {
        hit <- which(lower == tolower(cand))
        if (length(hit) > 0) return(nms[hit[1]])
    }

    NA_character_
}

clean_values <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == "." | x == "none" | x == "NA"] <- "NA"
    x
}

as_num <- function(x) {
    suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

write_count_table <- function(path, counts) {
    df <- data.frame(
        category = names(counts),
        count = as.integer(counts),
        stringsAsFactors = FALSE
    )
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

count_vec <- function(x) {
    sort(table(clean_values(x)), decreasing = TRUE)
}

top_counts <- function(counts, n = 25) {
    counts <- sort(counts, decreasing = TRUE)
    if (length(counts) <= n) return(counts)
    head(counts, n)
}

plot_horizontal_counts <- function(counts, main, xlab = "Count", top_n = NULL) {
    if (!is.null(top_n)) counts <- top_counts(counts, top_n)

    counts <- counts[order(counts, decreasing = FALSE)]

    par(mar = c(5, 15, 4, 3))
    bp <- barplot(
        counts,
        horiz = TRUE,
        las = 1,
        main = main,
        xlab = xlab,
        cex.names = 0.72
    )

    text(
        x = as.numeric(counts),
        y = bp,
        labels = as.integer(counts),
        pos = 4,
        cex = 0.7,
        xpd = TRUE
    )
}

make_plot <- function(fig, title, stub, plot_fun, width = 2400, height = 1800, res = 220) {
    png_path <- file.path(FIGURE_DIR, paste0(stub, ".png"))
    pdf_path <- file.path(FIGURE_DIR, paste0(stub, ".pdf"))

    status <- "created"
    notes <- ""

    tryCatch({
        png(png_path, width = width, height = height, res = res)
        plot_fun()
        dev.off()

        pdf(pdf_path, width = width / res, height = height / res)
        plot_fun()
        dev.off()
    }, error = function(e) {
        status <<- "failed"
        notes <<- conditionMessage(e)
        while (dev.cur() > 1) dev.off()
    })

    add_manifest(fig, title, status, png_path, pdf_path, notes)
}

skip_plot <- function(fig, title, notes) {
    message("Skipping ", fig, ": ", notes)
    add_manifest(fig, title, "skipped", "", "", notes)
}

parse_busco_status <- function(path) {
    if (is.na(path) || path == "" || path == "NA" || !file.exists(path)) {
        return(NULL)
    }

    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    lines <- lines[!grepl("^#", lines)]

    if (length(lines) == 0) return(NULL)

    split_lines <- strsplit(lines, "\t", fixed = TRUE)
    split_lines <- split_lines[lengths(split_lines) >= 2]

    if (length(split_lines) == 0) return(NULL)

    busco_id <- vapply(split_lines, function(x) x[1], character(1))
    status_raw <- vapply(split_lines, function(x) x[2], character(1))

    normalise_status <- function(x) {
        x <- tolower(x)
        if (grepl("duplicated", x)) return("Complete duplicated")
        if (grepl("complete", x)) return("Complete single-copy")
        if (grepl("fragmented", x)) return("Fragmented")
        if (grepl("missing", x)) return("Missing")
        "Other"
    }

    status <- vapply(status_raw, normalise_status, character(1))
    by_id <- split(status, busco_id)

    final_status <- vapply(by_id, function(s) {
        if ("Complete duplicated" %in% s) return("Complete duplicated")
        if ("Complete single-copy" %in% s) return("Complete single-copy")
        if ("Fragmented" %in% s) return("Fragmented")
        if ("Missing" %in% s) return("Missing")
        s[1]
    }, character(1))

    sort(table(final_status), decreasing = TRUE)
}

bin_exons <- function(x) {
    x <- as_num(x)
    x <- x[is.finite(x)]

    cut(
        x,
        breaks = c(-Inf, 1, 5, 10, Inf),
        labels = c("single-exon", "2-5 exons", "6-10 exons", ">10 exons"),
        right = TRUE
    )
}

bin_cluster_size <- function(x) {
    x <- as_num(x)
    x <- x[is.finite(x)]

    cut(
        x,
        breaks = c(-Inf, 1, 2, 4, 9, Inf),
        labels = c("1", "2", "3-4", "5-9", "10+"),
        right = TRUE
    )
}

master <- read_tsv(MASTER_TSV, "master annotation TSV")
dark <- read_tsv(DARK_TSV, "dark candidate TSV")
linked <- read_tsv(GENOME_LINKED_TSV, "genome-linked candidate TSV")
scaffold_summary <- read_tsv(SCAFFOLD_SUMMARY_TSV, "scaffold summary TSV")
busco_scaffold <- read_tsv(BUSCO_SCAFFOLD_TSV, "BUSCO scaffold summary TSV")
clusters <- read_tsv(CLUSTER_TSV, "dark candidate protein clusters TSV")
prioritised <- read_tsv(PRIORITISED_TSV, "prioritised candidate TSV")

# Figure 02 -------------------------------------------------------------------
fig <- "Figure 02"
title <- "Overall annotation-class breakdown across all proteins"

if (!is.null(master)) {
    class_col <- detect_col(
        master,
        c("annotation_class", "dark_candidate_category", "classification",
          "function_classification", "functional_classification", "class")
    )

    if (!is.na(class_col)) {
        counts <- count_vec(master[[class_col]])

        write_count_table(
            file.path(TABLE_DIR, "figure_02_annotation_class_counts.tsv"),
            counts
        )

        make_plot(
            fig,
            title,
            "figure_02_annotation_class_breakdown",
            function() {
                plot_horizontal_counts(
                    counts,
                    main = title,
                    xlab = "Number of proteins",
                    top_n = 30
                )
            }
        )
    } else {
        skip_plot(fig, title, "No annotation/classification column detected in master TSV.")
    }
} else {
    skip_plot(fig, title, "Master annotation TSV missing or unreadable.")
}

# Figure 03 -------------------------------------------------------------------
fig <- "Figure 03"
title <- "Dark candidate class composition"

if (!is.null(dark)) {
    class_col <- detect_col(
        dark,
        c("dark_candidate_category", "annotation_class", "classification",
          "function_classification", "functional_classification", "class")
    )

    if (!is.na(class_col)) {
        counts <- table(clean_values(dark[[class_col]]))

        target_classes <- c(
            "function_dark_no_current_annotation",
            "function_dark_but_signalp_secretory_candidate"
        )

        for (target in target_classes) {
            if (!(target %in% names(counts))) {
                counts <- c(counts, setNames(0, target))
            }
        }

        counts <- sort(counts, decreasing = TRUE)

        write_count_table(
            file.path(TABLE_DIR, "figure_03_dark_candidate_class_counts.tsv"),
            counts
        )

        make_plot(
            fig,
            title,
            "figure_03_dark_candidate_class_composition",
            function() {
                plot_horizontal_counts(
                    counts,
                    main = title,
                    xlab = "Number of dark candidates"
                )
            }
        )
    } else {
        skip_plot(fig, title, "No candidate class column detected.")
    }
} else {
    skip_plot(fig, title, "Dark candidate TSV missing or unreadable.")
}

# Figure 04 -------------------------------------------------------------------
fig <- "Figure 04"
title <- "Dark candidate FASTA sequence recovery"

seq_df <- NULL

if (!is.null(dark) && "sequence_status" %in% names(dark)) {
    seq_df <- dark
} else if (!is.null(linked) && "sequence_status" %in% names(linked)) {
    seq_df <- linked
}

if (!is.null(seq_df)) {
    counts <- count_vec(seq_df[["sequence_status"]])

    write_count_table(
        file.path(TABLE_DIR, "figure_04_fasta_recovery_counts.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_04_fasta_recovery",
        function() {
            plot_horizontal_counts(
                counts,
                main = title,
                xlab = "Number of dark candidates"
            )
        }
    )
} else {
    skip_plot(fig, title, "No sequence_status column found in dark or linked TSV.")
}

# Figure 05 -------------------------------------------------------------------
fig <- "Figure 05"
title <- "Genome-linking success of dark candidates"

if (!is.null(linked) && "genome_match_status" %in% names(linked)) {
    counts <- count_vec(linked[["genome_match_status"]])

    write_count_table(
        file.path(TABLE_DIR, "figure_05_genome_linking_status_counts.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_05_genome_linking_success",
        function() {
            plot_horizontal_counts(
                counts,
                main = title,
                xlab = "Number of dark candidate records"
            )
        }
    )
} else {
    skip_plot(fig, title, "genome_match_status column not found.")
}

# Figure 06 -------------------------------------------------------------------
fig <- "Figure 06"
title <- "Dark candidates per scaffold"

scaf_counts <- NULL

if (!is.null(scaffold_summary) &&
    all(c("contig", "n_dark_candidate_ids") %in% names(scaffold_summary))) {

    scaf_counts <- setNames(
        as_num(scaffold_summary[["n_dark_candidate_ids"]]),
        scaffold_summary[["contig"]]
    )
    scaf_counts <- scaf_counts[is.finite(scaf_counts)]

} else if (!is.null(linked) && "contig" %in% names(linked)) {
    tmp <- linked[linked[["contig"]] != "NA" & linked[["contig"]] != "", , drop = FALSE]
    scaf_counts <- table(tmp[["contig"]])
}

if (!is.null(scaf_counts) && length(scaf_counts) > 0) {
    scaf_counts <- sort(scaf_counts, decreasing = TRUE)

    write_count_table(
        file.path(TABLE_DIR, "figure_06_dark_candidates_per_scaffold.tsv"),
        scaf_counts
    )

    make_plot(
        fig,
        title,
        "figure_06_dark_candidates_per_scaffold_top25",
        function() {
            plot_horizontal_counts(
                scaf_counts,
                main = "Top 25 scaffolds by dark-candidate count",
                xlab = "Number of dark candidates",
                top_n = 25
            )
        }
    )
} else {
    skip_plot(fig, title, "No scaffold/count information available.")
}

# Figure 07 -------------------------------------------------------------------
fig <- "Figure 07"
title <- "Exon count distribution of genome-linked dark candidates"

if (!is.null(linked) && "exon_count" %in% names(linked)) {
    exon_bins <- bin_exons(linked[["exon_count"]])
    counts <- table(exon_bins)
    desired <- c("single-exon", "2-5 exons", "6-10 exons", ">10 exons")

    for (d in desired) {
        if (!(d %in% names(counts))) {
            counts <- c(counts, setNames(0, d))
        }
    }

    counts <- counts[desired]

    write_count_table(
        file.path(TABLE_DIR, "figure_07_exon_count_bins.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_07_exon_count_distribution",
        function() {
            par(mar = c(7, 5, 4, 2))
            bp <- barplot(
                counts,
                las = 2,
                main = title,
                ylab = "Number of dark candidates"
            )
            text(bp, counts, labels = counts, pos = 3, cex = 0.8)
        }
    )
} else {
    skip_plot(fig, title, "exon_count column not found.")
}

# Figure 08 -------------------------------------------------------------------
fig <- "Figure 08"
title <- "CDS and transcript span distributions"

if (!is.null(linked) &&
    all(c("cds_span_bp", "transcript_span_bp") %in% names(linked))) {

    cds <- as_num(linked[["cds_span_bp"]])
    tx <- as_num(linked[["transcript_span_bp"]])

    cds <- cds[is.finite(cds) & cds > 0]
    tx <- tx[is.finite(tx) & tx > 0]

    both <- linked[
        is.finite(as_num(linked[["cds_span_bp"]])) &
        is.finite(as_num(linked[["transcript_span_bp"]])) &
        as_num(linked[["cds_span_bp"]]) > 0 &
        as_num(linked[["transcript_span_bp"]]) > 0,
        ,
        drop = FALSE
    ]

    if (length(cds) > 0 && length(tx) > 0 && nrow(both) > 0) {
        span_table <- data.frame(
            metric = c("cds_span_bp", "transcript_span_bp"),
            n = c(length(cds), length(tx)),
            median_bp = c(median(cds), median(tx)),
            min_bp = c(min(cds), min(tx)),
            max_bp = c(max(cds), max(tx)),
            stringsAsFactors = FALSE
        )

        write.table(
            span_table,
            file.path(TABLE_DIR, "figure_08_span_distribution_summary.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        make_plot(
            fig,
            title,
            "figure_08_cds_and_transcript_span_distributions",
            function() {
                par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

                hist(
                    log10(cds),
                    breaks = 40,
                    main = "CDS span",
                    xlab = "log10(CDS span, bp)"
                )

                hist(
                    log10(tx),
                    breaks = 40,
                    main = "Transcript span",
                    xlab = "log10(transcript span, bp)"
                )

                plot(
                    log10(as_num(both[["transcript_span_bp"]])),
                    log10(as_num(both[["cds_span_bp"]])),
                    pch = 16,
                    cex = 0.5,
                    xlab = "log10(transcript span, bp)",
                    ylab = "log10(CDS span, bp)",
                    main = "CDS vs transcript span"
                )
                abline(0, 1, lty = 2)
            },
            width = 3300,
            height = 1300,
            res = 220
        )
    } else {
        skip_plot(fig, title, "No positive finite CDS/transcript spans available.")
    }
} else {
    skip_plot(fig, title, "cds_span_bp and/or transcript_span_bp columns not found.")
}

# Figure 09 -------------------------------------------------------------------
fig <- "Figure 09"
title <- "BUSCO status summary"

busco_counts <- parse_busco_status(BUSCO_FULL_TABLE)

if (!is.null(busco_counts) && length(busco_counts) > 0) {
    desired <- c("Complete single-copy", "Complete duplicated", "Fragmented", "Missing", "Other")

    for (x in desired) {
        if (!(x %in% names(busco_counts))) {
            busco_counts <- c(busco_counts, setNames(0, x))
        }
    }

    busco_counts <- busco_counts[desired]

    write_count_table(
        file.path(TABLE_DIR, "figure_09_busco_status_counts.tsv"),
        busco_counts
    )

    make_plot(
        fig,
        title,
        "figure_09_busco_status_summary",
        function() {
            par(mar = c(7, 5, 4, 2))
            bp <- barplot(
                busco_counts,
                las = 2,
                main = title,
                ylab = "Number of BUSCO groups"
            )
            text(bp, busco_counts, labels = busco_counts, pos = 3, cex = 0.8)
        }
    )
} else {
    skip_plot(fig, title, "BUSCO full_table.tsv not found or could not be parsed.")
}

# Figure 10 -------------------------------------------------------------------
fig <- "Figure 10"
title <- "Duplicated BUSCO burden per scaffold"

if (!is.null(busco_scaffold) &&
    all(c("contig", "duplicated_BUSCO_unique_ids") %in% names(busco_scaffold))) {

    counts <- setNames(
        as_num(busco_scaffold[["duplicated_BUSCO_unique_ids"]]),
        busco_scaffold[["contig"]]
    )
    counts <- counts[is.finite(counts)]
    counts <- sort(counts, decreasing = TRUE)

    write_count_table(
        file.path(TABLE_DIR, "figure_10_duplicated_buscos_per_scaffold.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_10_duplicated_buscos_per_scaffold_top25",
        function() {
            plot_horizontal_counts(
                counts,
                main = "Top 25 scaffolds by duplicated BUSCO count",
                xlab = "Unique duplicated BUSCO IDs",
                top_n = 25
            )
        }
    )
} else {
    skip_plot(fig, title, "BUSCO scaffold duplication summary missing or incomplete.")
}

# Figure 11 -------------------------------------------------------------------
fig <- "Figure 11"
title <- "Dark candidates versus duplicated BUSCO burden per scaffold"

dark_scaf_df <- NULL

if (!is.null(scaffold_summary) &&
    all(c("contig", "n_dark_candidate_ids") %in% names(scaffold_summary))) {

    dark_scaf_df <- data.frame(
        contig = scaffold_summary[["contig"]],
        n_dark_candidate_ids = as_num(scaffold_summary[["n_dark_candidate_ids"]]),
        stringsAsFactors = FALSE
    )

} else if (!is.null(linked) && "contig" %in% names(linked)) {
    tmp <- table(linked[["contig"]])
    dark_scaf_df <- data.frame(
        contig = names(tmp),
        n_dark_candidate_ids = as.integer(tmp),
        stringsAsFactors = FALSE
    )
}

if (!is.null(dark_scaf_df) &&
    !is.null(busco_scaffold) &&
    all(c("contig", "duplicated_BUSCO_unique_ids") %in% names(busco_scaffold))) {

    busco_df <- data.frame(
        contig = busco_scaffold[["contig"]],
        duplicated_BUSCO_unique_ids = as_num(busco_scaffold[["duplicated_BUSCO_unique_ids"]]),
        stringsAsFactors = FALSE
    )

    merged <- merge(dark_scaf_df, busco_df, by = "contig", all.x = TRUE)
    merged[["duplicated_BUSCO_unique_ids"]][is.na(merged[["duplicated_BUSCO_unique_ids"]])] <- 0
    merged <- merged[
        is.finite(merged[["n_dark_candidate_ids"]]) &
        is.finite(merged[["duplicated_BUSCO_unique_ids"]]),
        ,
        drop = FALSE
    ]

    write.table(
        merged,
        file.path(TABLE_DIR, "figure_11_dark_candidates_vs_duplicated_buscos.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    make_plot(
        fig,
        title,
        "figure_11_dark_candidates_vs_duplicated_buscos",
        function() {
            x <- merged[["duplicated_BUSCO_unique_ids"]]
            y <- merged[["n_dark_candidate_ids"]]

            plot(
                x,
                y,
                pch = 16,
                cex = 0.8,
                xlab = "Unique duplicated BUSCO IDs per scaffold",
                ylab = "Dark candidates per scaffold",
                main = title
            )

            fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
            if (!is.null(fit)) abline(fit, lty = 2)

            top <- order(y + x, decreasing = TRUE)
            top <- head(top, min(10, length(top)))

            text(
                x[top],
                y[top],
                labels = substr(merged[["contig"]][top], 1, 22),
                pos = 4,
                cex = 0.65,
                xpd = TRUE
            )
        }
    )
} else {
    skip_plot(fig, title, "Need both scaffold dark-candidate counts and BUSCO scaffold summary.")
}

# Figure 12 -------------------------------------------------------------------
fig <- "Figure 12"
title <- "Near-identical dark-candidate protein cluster-size distribution"

if (!is.null(clusters) &&
    all(c("dark_cluster_id", "cluster_size") %in% names(clusters))) {

    cluster_unique <- clusters[!duplicated(clusters[["dark_cluster_id"]]), , drop = FALSE]
    bins <- bin_cluster_size(cluster_unique[["cluster_size"]])
    counts <- table(bins)

    desired <- c("1", "2", "3-4", "5-9", "10+")
    for (d in desired) {
        if (!(d %in% names(counts))) {
            counts <- c(counts, setNames(0, d))
        }
    }

    counts <- counts[desired]

    write_count_table(
        file.path(TABLE_DIR, "figure_12_cluster_size_bins.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_12_dark_candidate_cluster_size_distribution",
        function() {
            par(mar = c(6, 5, 4, 2))
            bp <- barplot(
                counts,
                las = 2,
                main = title,
                ylab = "Number of clusters"
            )
            text(bp, counts, labels = counts, pos = 3, cex = 0.8)
        }
    )
} else {
    skip_plot(fig, title, "Cluster table missing or required columns absent.")
}

# Figure 13 -------------------------------------------------------------------
fig <- "Figure 13"
title <- "Duplication interpretation categories"

if (!is.null(prioritised) && "duplication_interpretation" %in% names(prioritised)) {
    counts <- count_vec(prioritised[["duplication_interpretation"]])

    write_count_table(
        file.path(TABLE_DIR, "figure_13_duplication_interpretation_counts.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_13_duplication_interpretation_categories",
        function() {
            plot_horizontal_counts(
                counts,
                main = title,
                xlab = "Number of dark candidate records"
            )
        }
    )
} else {
    skip_plot(fig, title, "duplication_interpretation column not found in prioritised table.")
}

# Figure 14 -------------------------------------------------------------------
fig <- "Figure 14"
title <- "Final priority tier summary"

if (!is.null(prioritised) && "priority_tier" %in% names(prioritised)) {
    counts <- count_vec(prioritised[["priority_tier"]])

    preferred <- c("high_priority", "medium_priority", "low_priority_or_manual_review")
    for (x in preferred) {
        if (!(x %in% names(counts))) {
            counts <- c(counts, setNames(0, x))
        }
    }

    counts <- counts[preferred]

    write_count_table(
        file.path(TABLE_DIR, "figure_14_priority_tier_counts.tsv"),
        counts
    )

    make_plot(
        fig,
        title,
        "figure_14_priority_tier_summary",
        function() {
            par(mar = c(8, 5, 4, 2))
            bp <- barplot(
                counts,
                las = 2,
                main = title,
                ylab = "Number of dark candidate records"
            )
            text(bp, counts, labels = counts, pos = 3, cex = 0.8)
        }
    )
} else {
    skip_plot(fig, title, "priority_tier column not found in prioritised table.")
}

# Figure 15 -------------------------------------------------------------------
fig <- "Figure 15"
title <- "Priority score distribution"

if (!is.null(prioritised) && "priority_score" %in% names(prioritised)) {
    scores <- as_num(prioritised[["priority_score"]])
    scores <- scores[is.finite(scores)]

    if (length(scores) > 0) {
        score_table <- data.frame(
            priority_score = sort(unique(scores)),
            count = as.integer(table(scores)[as.character(sort(unique(scores)))]),
            stringsAsFactors = FALSE
        )

        write.table(
            score_table,
            file.path(TABLE_DIR, "figure_15_priority_score_counts.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        make_plot(
            fig,
            title,
            "figure_15_priority_score_distribution",
            function() {
                hist(
                    scores,
                    breaks = seq(min(scores) - 0.5, max(scores) + 0.5, by = 1),
                    main = title,
                    xlab = "Priority score",
                    ylab = "Number of dark candidate records"
                )
            }
        )
    } else {
        skip_plot(fig, title, "priority_score column contained no finite numeric values.")
    }
} else {
    skip_plot(fig, title, "priority_score column not found in prioritised table.")
}

manifest_path <- file.path(FIGURE_DIR, "figure_manifest.tsv")

write.table(
    manifest,
    manifest_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

message("Figure generation complete.")
message("Figure directory: ", FIGURE_DIR)
message("Figure tables: ", TABLE_DIR)
message("Manifest: ", manifest_path)
RSCRIPT

chmod +x "${RSCRIPT_PATH}"

echo
echo "Running R figure-generation script..."
echo "Project directory: ${PROJECT_DIR}"
echo "Figure directory: ${FIGURE_DIR}"
echo "Master TSV: ${MASTER_TSV}"
echo "Dark TSV: ${DARK_TSV}"
echo "Genome-linked TSV: ${GENOME_LINKED_TSV}"
echo "Scaffold summary TSV: ${SCAFFOLD_SUMMARY_TSV}"
echo "BUSCO full table: ${BUSCO_FULL_TABLE}"
echo "BUSCO scaffold TSV: ${BUSCO_SCAFFOLD_TSV}"
echo "Cluster TSV: ${CLUSTER_TSV}"
echo "Prioritised TSV: ${PRIORITISED_TSV}"
echo

Rscript "${RSCRIPT_PATH}"

echo
echo "Generated figure files:"
find "${FIGURE_DIR}" -maxdepth 1 -type f \( -name "*.png" -o -name "*.pdf" -o -name "figure_manifest.tsv" \) | sort

echo
echo "Manifest:"
cat "${FIGURE_DIR}/figure_manifest.tsv"