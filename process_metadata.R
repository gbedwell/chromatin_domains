library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

dat <- data.table::fread("son_cleaned.tsv", header = FALSE)

links <- dat$V1

input <- dat |>
    dplyr::filter(grepl("input", dat$V24) | grepl("gDNA", dat$V24)) |>
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
    dplyr::filter(!grepl("input", dat$V24) & !grepl("gDNA", dat$V24)) |>
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
            file="son_tsa_download_links.txt",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

write.table(dat.matched,
            file="son_tsa_file_match.txt",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

rm()

dat <- data.table::fread("laminB1_cleaned.tsv", header = FALSE)

links <- dat$V1

input <- dat |>
    dplyr::filter(grepl("input", dat$V24) | grepl("gDNA", dat$V24)) |>
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
    dplyr::filter(!grepl("input", dat$V24) & !grepl("gDNA", dat$V24)) |>
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
            file="laminB1_tsa_download_links.txt",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

write.table(dat.matched,
            file="laminB1_tsa_file_match.txt",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)
