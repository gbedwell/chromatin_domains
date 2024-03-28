#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

dat <- data.table::fread(paste0( args[1], "_cleaned.tsv" ), header = FALSE)

links <- dat$V1

input <- dat |>
    dplyr::filter(grepl("input", dat$V24) | grepl("gDNA", dat$V24) | grepl("None", dat$V21)) |>
    dplyr::arrange(V22) |>
    dplyr::mutate(V22 = gsub(" .*", "", V22)) |>
    dplyr::select(V1, V2, V3, V4, V9, V22, V24) |>
    dplyr::mutate(
        V9 = case_when(
            grepl("shock", V24) ~ 2,
            TRUE ~ V9
        ),
        V22 = paste0(V22, "_input_replicate_", V9)
    ) |>
    dplyr::select(V4, V22)

pulldown <- dat |>
    dplyr::filter(!grepl("input", dat$V24) & !grepl("gDNA", dat$V24) & !grepl("None", dat$V21)) |>
    dplyr::arrange(V22) |>
    dplyr::mutate(V22 = gsub(" .*", "", V22)) |>
    dplyr::select(V1, V2, V3, V4, V9, V22, V24) |>
    dplyr::mutate(
        V9 = case_when(
            grepl("shock", V24) ~ 2,
            TRUE ~ V9
        ),
        V22 = paste0(V22, "_pulldown_replicate_", V9)
    ) |>
    dplyr::select(V4, V22)

dat.matched <- rbind(input, pulldown) |>
    dplyr::arrange(V22)

write.table(links,
            file=paste0( args[1], "_tsa_download_links.txt" ),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

write.table(dat.matched,
            file=paste0( args[1], "_tsa_file_match.txt" ),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

rm()