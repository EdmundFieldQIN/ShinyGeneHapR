#' @name HaploWorldMap
#' @title Display Haplotype Distribution on World Map using ggplot2
#' @description
#' Create a world map visualization showing the geographic distribution of
#' haplotypes using ggplot2 and scatterpie. This function generates a
#' publication-quality map with pie charts representing haplotype frequencies
#' at different geographic locations.
#'
#' @importFrom ggplot2 ggplot geom_polygon aes coord_quickmap theme_bw
#' @importFrom ggplot2 element_rect element_blank element_text annotate
#' @importFrom ggplot2 theme scale_fill_manual
#' @importFrom dplyr filter mutate
#' @importFrom scatterpie geom_scatterpie geom_scatterpie_legend
#' @importFrom reshape2 dcast
#' @importFrom utils askYesNo
#'
#' @examples
#' \donttest{
#' data("geneHapR_test")
#' HaploWorldMap(
#'    hap = hapResult,
#'    AccINFO = AccINFO,
#'    LON.col = "Longitude",
#'    LAT.col = "Latitude",
#'    hapNames = c("H001", "H002", "H003"),
#'    hapColor = c("#F2406D", "#0084FF", "#00E3BA"),
#'    pie.log.base = 10,
#'    pie.line.color = "white",
#'    pie.line.width = 0,
#'    pie.legend.size = 2,
#'    pie.legend.posi = c(-150, -30),
#'    pie.legend.title.size = 5,
#'    pie.legend.title.posi = c(-140, -20),
#'    hap.legend.size = 10,
#'    hap.legend.posi = "right",
#'    map.fill.color = "grey60",
#'    map.border.color = "white",
#'    map.border.width = 0.2,
#'    map.background.color = "white",
#'    exclude.regions = TRUE,
#'    exclude.regions.names = "Antarctica"
#' )
#'
#' @param hap an object of hapResult or hapSummary class containing
#' haplotype information
#' @param AccINFO a data.frame containing accession information including
#' longitude and latitude coordinates
#' @param LON.col column name for longitude coordinates in AccINFO
#' @param LAT.col column name for latitude coordinates in AccINFO
#' @param hapNames vector of haplotype names to display on the map
#' @param hapColor optional vector of colors for each haplotype. If NULL,
#' default ggplot2 colors will be used
#' @param pie.log.base base for logarithmic scaling of pie chart sizes
#' (default: 10)
#' @param pie.line.color color of pie chart borders
#' @param pie.line.width line width of pie chart borders (0 for no border)
#' @param pie.legend.size size of the pie size legend
#' @param pie.legend.posi x,y coordinates for the pie size legend position
#' @param pie.legend.title.size text size for the pie size legend title
#' @param pie.legend.title.posi x,y coordinates for the pie size legend title
#' @param hap.legend.size text size for the haplotype legend
#' @param hap.legend.posi position of the haplotype legend ("right", "left",
#' "top", "bottom", or a numeric vector of x,y coordinates)
#' @param map.fill.color fill color for map polygons
#' @param map.border.color border color for map polygons
#' @param map.border.width border line width for map polygons
#' @param map.background.color background color for the map panel
#' @param exclude.regions logical, whether to exclude specific regions from
#' the world map
#' @param exclude.regions.names names of regions to exclude (default: "Antarctica")
#' @param ... additional arguments passed to ggplot2 theme elements
#'
#' @return A ggplot2 object containing the haplotype world map visualization.
#' The plot can be further customized using standard ggplot2 functions.
#'
#' @export
#'
HaploWorldMap <- function(
    hap,
    AccINFO,
    LON.col = "Longitude",
    LAT.col = "Latitude",
    hapNames,
    hapColor = NULL,
    pie.log.base = 10,
    pie.line.color = "black",
    pie.line.width = 0,
    pie.legend.size = 2,
    pie.legend.posi = c(-150, -30),
    pie.legend.title.size = 5,
    pie.legend.title.posi = c(-140, -20),
    hap.legend.size = 10,
    hap.legend.posi = "right",
    map.fill.color = "grey60",
    map.border.color = "white",
    map.border.width = 0.2,
    map.background.color = "white",
    exclude.regions = TRUE,
    exclude.regions.names = "Antarctica",
    ...
) {
    hap <- na.omit(hap)
    # 判断hap对象类型并提取样本到haplotype的映射关系
    if (inherits(hap, 'hapResult')) {
        # 如果是hapResult类型：直接提取Accession和Hap列
        acc2hap <- hap[, "Hap"]
        names(acc2hap) <- hap[, "Accession"]
    } else if (inherits(hap, "hapSummary")) {
        # 如果是hapSummary类型：需要根据频率展开
        acc2hap <- rep(hap[, "Hap"], hap[, "freq"])
        names(acc2hap) <- unlist(strsplit(hap[, "Accession"], ";"))
    }
    # 将haplotype信息添加到AccINFO中
    AccINFO$Hap <- acc2hap[row.names(AccINFO)]
    # 提取地理数据：Hap、经度、纬度
    geoData <- AccINFO[c("Hap", LON.col, LAT.col)]
    geoData$value <- 1 # 添加一个值列用于计数
    # Check the column data format
    if (!inherits(geoData[, LAT.col], "numeric")) {
        warning("The '", LAT.col, "' column should be numeric")
        geoData[, LAT.col] <- suppressWarnings(as.numeric(geoData[, LAT.col]))
        if (!inherits(geoData[, LAT.col], "numeric")) {
            stop("We couldn't conver the '", LAT.col, "' column as 'numeric'")
        }
    }
    if (!inherits(geoData[, LON.col], "numeric")) {
        warning("The '", LON.col, "' column should be numeric")
        geoData[, LON.col] <- suppressWarnings(as.numeric(geoData[, LON.col]))
        if (!inherits(geoData[, LON.col], "numeric")) {
            stop("We couldn't conver the '", LON.col, "' column as 'numeric'")
        }
    }
    # 使用reshape2::dcast将数据转换为宽格式
    formu <- paste0(LON.col, "+", LAT.col, "~Hap")
    dF <- reshape2::dcast(
        geoData,
        formula = formu,
        value.var = "value",
        fun.aggregate = length
    )

    # 检查所有指定的haplotype是否在数据中
    if (!all(hapNames %in% names(dF))) {
        # 如果某些haplotype不存在，询问用户是否继续
        n_names <- hapNames[!hapNames %in% names(dF)]
        proceed <- utils::askYesNo(
            paste0(
                "We could not find any individual belong to '",
                n_names,
                "' in your AccINFO, Procceed?"
            ),
            default = FALSE
        )
        if (proceed) {
            hapNames <- hapNames[hapNames %in% names(dF)]
        } else {
            return(NULL)
        }
    }

    # 筛选dF中只包含hapNames列的数据，生成新的数据，并在最后一列加上Size列，Size为当前行hapNames列的和
    dF <- dF[, c(LON.col, LAT.col, hapNames)] %>%
        dplyr::mutate(Size = rowSums(.[, hapNames]))
    dF <- na.omit(dF)
    # 读取世界地图数据------
    world_map <- ggplot2::map_data("world")
    # 判断是否需要排除特定区域
    if (exclude.regions) {
        world_map <- dplyr::filter(world_map, !region %in% exclude.regions.names)
    }

    p <- ggplot() +
        ggplot2::geom_polygon(
            data = world_map,
            aes(x = long, y = lat, group = group),
            color = map.border.color, # 地图边界颜色
            fill = map.fill.color, # 地图填充颜色
            linewidth = map.border.width # 地图边界粗细
        ) +
        coord_quickmap() +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = map.background.color),
            panel.grid = element_blank()
        ) +
        scatterpie::geom_scatterpie(
            data = dF,
            aes(x = Longitude, y = Latitude, r = log(base = pie.log.base, x = Size + 1) * 2),
            cols = colnames(dF[, 3:(length(hapNames) + 2)]),
            color = pie.line.color,
            linewidth = pie.line.width,
            legend_name = "Haplotype"
        ) +
        scatterpie::geom_scatterpie_legend(
            r = log(base = pie.log.base, x = dF$Size + 1) * 2,
            x = pie.legend.posi[1],
            y = pie.legend.posi[2],
            # 将图例标签转换为原始的Size值（即log(base = 10, x + 1) * 2的反函数）
            labeller = function(r) round(pie.log.base^(r / 2) - 1, 0),
            size = pie.legend.size
        ) +
        annotate("text", x = pie.legend.title.posi[1], y = pie.legend.title.posi[2], label = "Accessions", size = pie.legend.title.size) +
        theme(
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.text = element_text(size = hap.legend.size, colour = 'grey20'),
            legend.position = hap.legend.posi
        )
    if (!is.null(hapColor)) {
        p <- p + scale_fill_manual(values = hapColor)
    }
    return(p)
}
