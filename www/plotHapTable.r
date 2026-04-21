plotHapTable3 <- function(
    hapSummary,
    hapPrefix = "H",
    title = "",
    geneName = geneName,
    INFO_tag = NULL,
    tag_split = tag_split,
    tag_field = tag_field,
    tag_name = tag_name,
    displayIndelSize = 0,
    angle = c(0, 45, 90),
    replaceMultiAllele = TRUE,
    ALLELE.color = "grey90",
    genotype_colors = NULL, # 命名向量，如 c("A"="#1b9e77","T"="#d95f02","C"="#7570b3","G"="#e7298a")
    text_size = 3.5,
    axis_text_size = 11,
    title_size = 14,
    caption_size = 9,
    font_family = "sans",
    cell_border_color = "white",
    cell_border_size = 0.4,
    show_legend = TRUE,
    legend_title = "Genotype"
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package 'tidyr' is required.")
    }
    if (!requireNamespace("stringr", quietly = TRUE)) {
        stop("Package 'stringr' is required.")
    }
    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Package 'scales' is required.")
    }

    # ---- 输入对象兼容 ----
    if (!inherits(hapSummary, "hapSummary")) {
        if (inherits(hapSummary, "hapResult")) {
            hapSummary <- hap_summary(hapSummary)
        } else {
            stop("'hapSummary' must inherit from class 'hapSummary' or 'hapResult'.")
        }
    }

    hapSummary <- as.data.frame(hapSummary, stringsAsFactors = FALSE)

    if (!missing(geneName) && identical(title, "")) {
        title <- geneName
    }

    if ("Accession" %in% colnames(hapSummary)) {
        hapSummary <- hapSummary[, colnames(hapSummary) != "Accession", drop = FALSE]
    }

    row_id_col <- colnames(hapSummary)[1]
    if (!"ALLELE" %in% hapSummary[[row_id_col]]) {
        stop("No 'ALLELE' row found in hapSummary.")
    }

    # ---- 提取hap行 ----
    hap_idx <- grepl(paste0("^", hapPrefix, "[0-9]+$"), hapSummary[[row_id_col]])
    hps <- as.matrix(hapSummary[hap_idx, , drop = FALSE])
    if (nrow(hps) < 1) {
        stop("No haplotype rows found. Please check 'hapSummary' and 'hapPrefix'.")
    }

    ALLELE <- unlist(hapSummary[hapSummary[[row_id_col]] == "ALLELE", , drop = TRUE], use.names = FALSE)

    # ---- Indel处理 ----
    footi <- ""
    indel_map <- character(0)

    probe_indel <- is.indel.allele(ALLELE)
    if (TRUE %in% probe_indel) {
        display_threshold <- displayIndelSize + 1
        all_indel <- ALLELE[probe_indel]
        all_indel <- unlist(stringr::str_split(all_indel, "[,/]"))
        all_indel <- unique(all_indel[stringr::str_length(all_indel) > display_threshold])

        if (length(all_indel) > 0 && !is.na(all_indel[1])) {
            indel_notes <- paste0("i", seq_len(length(all_indel)))
            names(indel_notes) <- all_indel
            indel_map <- indel_notes

            # 替换hap矩阵里的indel
            hps[hps %in% all_indel] <- indel_notes[hps[hps %in% all_indel]]

            # 替换ALLELE行里的indel
            for (i in seq_along(ALLELE)) {
                if (probe_indel[i]) {
                    token <- unlist(stringr::str_split(ALLELE[i], "[,/]"))
                    idx <- token %in% names(indel_notes)
                    token[idx] <- indel_notes[token[idx]]
                    ALLELE[i] <- gsub(",", "/", paste(token, collapse = ","), fixed = TRUE)
                }
            }

            footi <- paste(indel_notes, names(indel_notes), sep = ":", collapse = "; ")
        }
    }

    # ---- 多等位替换 ----
    footT <- ""
    multiallele_map <- character(0)
    probe_mula <- is.multiallelic.allele(ALLELE)
    if (isTRUE(replaceMultiAllele) && (TRUE %in% probe_mula)) {
        rept <- ALLELE[probe_mula]
        noteT <- paste0("T", seq_len(length(rept)))
        names(noteT) <- rept
        ALLELE[probe_mula] <- noteT
        multiallele_map <- noteT

        if (!is.na(names(noteT)[1])) {
            footT <- paste(noteT, names(noteT), sep = ":", collapse = "; ")
        }
    }

    foot <- if (nzchar(footi) && nzchar(footT)) {
        paste(footT, footi, sep = "\n")
    } else {
        paste0(footT, footi)
    }

    # 先插入ALLELE行
    hps <- rbind(ALLELE, hps)
    hps <- as.data.frame(hps, stringsAsFactors = FALSE)

    # ---- INFO_tag处理 ----
    if (!is.null(INFO_tag)) {
        m <- "length of 'tag_split', 'tag_name' and 'tag_field' should match 'INFO_tag'"

        if (!missing(tag_split) && length(INFO_tag) != length(tag_split)) {
            stop(m)
        }

        if (missing(tag_name)) {
            tag_name <- INFO_tag
        } else if (length(INFO_tag) != length(tag_name)) {
            stop(m)
        }

        info_rows <- hapSummary[hapSummary[[row_id_col]] == "INFO", , drop = FALSE]
        if (nrow(info_rows) == 0) {
            stop("'INFO_tag' provided but no INFO row found in hapSummary.")
        }

        INFO <- t(info_rows)
        INFO <- strsplit(INFO, ";")

        for (ntag in seq_along(INFO_tag)) {
            tag_i <- INFO_tag[ntag]

            INFOi <- unlist(lapply(INFO, function(x) {
                hit <- x[startsWith(x, paste0(tag_i, "="))]
                if (length(hit) == 0) NA_character_ else hit[1]
            }))
            INFOi <- stringr::str_remove(INFOi, paste0("^", tag_i, "="))

            if (missing(tag_field)) {
                Ii <- INFOi
            } else {
                if (tag_i %in% c("ANN", "SNPEFF")) {
                    if (missing(geneName)) {
                        stop("'geneName' is missing for ANN/SNPEFF parsing.")
                    }
                    Ii <- lapply(INFOi, function(x) {
                        if (is.na(x)) {
                            return(NA_character_)
                        }
                        ann_records <- strsplit(x, ",")[[1]]
                        ann_fields <- lapply(ann_records, function(y) strsplit(y, "[|]")[[1]])
                        vals <- vapply(
                            ann_fields,
                            function(z) {
                                if (length(z) >= 4 && identical(z[4], geneName) && length(z) >= tag_field[ntag]) {
                                    z[tag_field[ntag]]
                                } else {
                                    NA_character_
                                }
                            },
                            character(1)
                        )
                        vals <- vals[!is.na(vals) & nzchar(vals)]
                        if (length(vals) == 0) NA_character_ else paste(vals, collapse = ",")
                    })
                    Ii <- unlist(Ii)
                } else {
                    if (missing(tag_split)) {
                        stop("'tag_split' is missing.")
                    }
                    Ii <- unlist(lapply(INFOi, function(x) {
                        if (is.na(x)) {
                            return(NA_character_)
                        }
                        xs <- unlist(strsplit(x, tag_split))
                        if (length(xs) < tag_field[ntag]) NA_character_ else xs[tag_field[ntag]]
                    }))
                }
            }

            # 构建INFO行并按列名对齐
            info_vec <- rep(NA_character_, ncol(hps))
            names(info_vec) <- colnames(hps)
            info_vec[1] <- tag_name[ntag]

            fill_cols <- setdiff(colnames(hps), colnames(hps)[1])
            nfill <- min(length(fill_cols), length(Ii))
            if (nfill > 0) {
                info_vec[fill_cols[seq_len(nfill)]] <- Ii[seq_len(nfill)]
            }
            if ("freq" %in% colnames(hps)) {
                info_vec["freq"] <- ""
            }

            hps <- rbind(info_vec, hps)
            hps <- as.data.frame(hps, stringsAsFactors = FALSE)
        }
    }

    # ---- 转长表 ----
    colnames(hps)[1] <- "RowName"
    hps$RowName <- as.character(hps$RowName)

    long_df <- tidyr::pivot_longer(
        hps,
        cols = -RowName,
        names_to = "ColName",
        values_to = "Label"
    )

    long_df$Fill <- long_df$Label

    # 非基因型单元格设为NA，统一用ALLELE.color底色
    if (!missing(INFO_tag) && !is.null(INFO_tag)) {
        long_df$Fill[long_df$RowName %in% tag_name] <- NA_character_
    }
    long_df$Fill[long_df$RowName == "ALLELE"] <- NA_character_
    long_df$Fill[long_df$ColName == "freq"] <- NA_character_

    # 行顺序：Hn到H1，然后其他行
    row_levels_raw <- unique(long_df$RowName)
    hap_levels <- row_levels_raw[grepl(paste0("^", hapPrefix, "[0-9]+$"), row_levels_raw)]
    other_levels <- row_levels_raw[!row_levels_raw %in% hap_levels]
    row_levels <- c(hap_levels[order(hap_levels, decreasing = TRUE)], other_levels)
    long_df$RowName <- factor(long_df$RowName, levels = row_levels)

    # 角度参数
    angle <- angle[1]
    if (!angle %in% c(0, 45, 90)) {
        stop("angle should be one of 0, 45, and 90")
    }
    vjust <- if (angle == 0) {
        0.5
    } else if (angle == 45) {
        0.1
    } else {
        1
    }
    hjust <- if (angle == 0) {
        0.5
    } else if (angle == 45) {
        0.1
    } else {
        1
    }

    # ---- 颜色映射 ----
    genotype_levels <- sort(unique(stats::na.omit(long_df$Fill)))
    if (length(genotype_levels) > 0) {
        if (is.null(genotype_colors)) {
            pal <- scales::hue_pal()(length(genotype_levels))
            names(pal) <- genotype_levels
            genotype_colors <- pal
        } else {
            if (is.null(names(genotype_colors))) {
                if (length(genotype_colors) < length(genotype_levels)) {
                    stop("'genotype_colors' must be a named vector or have enough colors.")
                }
                names(genotype_colors) <- genotype_levels
                genotype_colors <- genotype_colors[genotype_levels]
            } else {
                miss <- setdiff(genotype_levels, names(genotype_colors))
                if (length(miss) > 0) {
                    extra <- scales::hue_pal()(length(miss))
                    names(extra) <- miss
                    genotype_colors <- c(genotype_colors, extra)
                }
                genotype_colors <- genotype_colors[genotype_levels]
            }
        }
    } else {
        genotype_colors <- c()
    }

    # 避免geom_text因NA warning
    long_df$LabelPlot <- ifelse(is.na(long_df$Label), "", as.character(long_df$Label))

    p <- ggplot2::ggplot(long_df, ggplot2::aes(x = ColName, y = RowName, fill = Fill)) +
        ggplot2::geom_tile(color = cell_border_color, linewidth = cell_border_size) +
        ggplot2::geom_text(
            ggplot2::aes(label = LabelPlot),
            size = text_size,
            family = font_family,
            na.rm = TRUE
        ) +
        ggplot2::scale_fill_manual(
            values = genotype_colors,
            na.value = ALLELE.color,
            drop = FALSE,
            guide = if (show_legend) "legend" else "none"
        ) +
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(position = "top")) +
        ggplot2::labs(
            title = title,
            caption = foot,
            fill = legend_title
        ) +
        ggplot2::theme_minimal(base_family = font_family) +
        ggplot2::theme(
            legend.position = if (show_legend) "right" else "none",
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(
                angle = angle,
                vjust = vjust,
                hjust = hjust,
                size = axis_text_size
            ),
            axis.text.y = ggplot2::element_text(size = axis_text_size),
            plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
            plot.caption = ggplot2::element_text(size = caption_size, hjust = 0)
        )

    attr(p, "indel_map") <- indel_map
    attr(p, "multiallele_map") <- multiallele_map
    p
}

# # 测试数据
# setwd("D:/R/R_Workspace/workspace/Haplotype_analysis/完全版")
# source(file.path(find.package("geneHapR"), "extdata", "www", "zzz.R"))
# AccINFO <- geneHapR::import_AccINFO(
#     "accinfo.csv",
#     sep = ",", # 分隔符号，默认为制表符"/t"
#     na.strings = "NA"
# )
# vcf <- geneHapR::import_vcf("LOC_Os05g46480.vcf")
# hapResult <- geneHapR::vcf2hap(vcf, hetero_remove = T, na_drop = T, hapPrefix = "H")
# hapSummary <- geneHapR::hap_summary(hapResult)

# p <- plotHapTable3(
#     hapSummary,
#     hapPrefix = "H",
#     title = "LOC_Os05g46480 Haplotype Table",
#     geneName = "LOC_Os05g46480",
#     INFO_tag = NULL,
#     tag_split = ",",
#     tag_field = 2,
#     tag_name = "Annotation",
#     displayIndelSize = 0,
#     angle = 45,
#     replaceMultiAllele = TRUE,
#     ALLELE.color = "#FFFFFF00",
#     genotype_colors = c("A" = "#F2406D", "T" = "#0084FF", "C" = "#7570b300", "G" = "#e7298b00"), # 命名向量，如 c("A"="#1b9e77","T"="#d95f02","C"="#7570b3","G"="#e7298a")
#     text_size = 3.5,
#     axis_text_size = 11,
#     title_size = 14,
#     caption_size = 9,
#     font_family = "sans",
#     cell_border_color = "white",
#     cell_border_size = 0.4,
#     show_legend = TRUE,
#     legend_title = "Genotype"
# )
# p
