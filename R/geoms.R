#' Plotting functions
#'
#' @description
#' 
#' Plotting functions for tidyCoverage objects. 
#' Two geoms are provided:
#' 
#' - `geom_coverage()`: for plotting coverages over individual loci.
#' - `geom_aggrcoverage()`: for plotting aggregated coverages with confidence intervals.
#' 
#' See the Details section for more information on the aesthetics used by each geom.
#' 
#' @details 
#' 
#' These geoms are drawn using `geom_line`/`ribbon`/`area()` so they support the same
#' aesthetics: `colour`, `linetype` and `linewidth.` Both geoms also support 
#' the `unit` argument to control the x axis units (b, kb, Mb). 
#' 
#' In addition, they each support additional arguments:
#' 
#' - `geom_coverage` uses a `type` argument to switch between line plot and area plots; 
#' - `geom_aggrcoverage` uses a `ci` argument to toggle the confidence interval display.
#' 
#' @name ggplot-tidyCoverage
#' @rdname ggplot-tidyCoverage
#' 
#' @param mapping Set of aesthetic mappings created by aes(). 
#'     By default, no color/fill aesthetic 
#'     is specified, but they can be assigned to a variable with `mapping = aes(...)`. 
#'     Note that `x` and `y` are automatically filled. 
#' @param data Data frame passed to geom_*. Typically a `CoverageExperiment` object 
#'     (expanded to a tibble) or a `AggregatedCoverage` object. 
#' @param type Choose between "line" and "area" style for `geom_coverage()`.
#' @param ci Should the confidence interval be plotted by `geom_aggrcoverage()`?
#'     (default: TRUE)
#' @param unit Rounding of x axis (any of c('b', 'kb', 'Mb')).
#' @param grid Should the plot grid by displayed? (default: FALSE).
#' @param alpha Transparency level for `geom_coverage()` (default: 0.6).
#' @param linewidth Line width for `geom_coverage()` (default: 0.4).
#' @param raster Should the plot be rasterized for faster rendering?
#'     (default: TRUE)
#' @param ...,na.rm,show.legend,inherit.aes Argument passed to `ggplot` 
#'     internal functions
#' @return A `ggplot` object
#'
#' @import ggplot2
#' @importFrom scales oob_squish
#' @importFrom scales unit_format
#'
#' @examples
#' library(rtracklayer)
#' library(plyranges)
#' library(ggplot2)
#' library(purrr)
#' TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
#' features <- list(
#'     TSS_fwd = import(TSSs_bed) |> filter(strand == '+'), 
#'     TSS_rev = import(TSSs_bed) |> filter(strand == '-'), 
#'     conv_sites = import(system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage"))
#' )
#' tracks <- list(
#'     RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
#'     RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage"), 
#'     Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage")
#' ) |> map(import, as = 'Rle')
#' ce <- CoverageExperiment(tracks, features, width = 5000, center = TRUE, scale = TRUE)
#' ac <- aggregate(ce) 
#' 
#' #############################################################################
#' ## 1. Plotting aggregated coverage
#' #############################################################################
#' 
#' ac |> 
#'     as_tibble() |> 
#'     ggplot() + 
#'     geom_aggrcoverage(aes(col = track)) + 
#'     facet_grid(track ~ features) + 
#'     geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', linewidth = 0.5)
#' 
#' #############################################################################
#' ## 2. Plotting track coverages over individual loci
#' #############################################################################
#' 
#' ce2 <- CoverageExperiment(
#'     tracks, 
#'     GRangesList(list(locus1 = "II:400001-455000", locus2 = "IV:720001-775000")), 
#'     window = 50
#' )
#' expand(ce2) |> 
#'     mutate(coverage = ifelse(track != 'Scc1', scales::oob_squish(coverage, c(0, 50)), coverage)) |>
#'     ggplot() + 
#'     geom_coverage(aes(fill = track)) + 
#'     facet_grid(track~features, scales = 'free')
NULL 

GeomAggrCoverage <- ggplot2::ggproto("GeomAggrCoverage", ggplot2::Geom,
    setup_params = function(data, params) {
        params$ci <- params$ci
        params
    },
    extra_params = c("na.rm"),
    required_aes = c("x", "y", "ymin", "ymax"), 
    default_aes = ggplot2::aes(
        colour = "black", 
        linewidth = 1, 
        linetype = 1, 
        alpha = 0.4
    ),
    draw_group = function(data, params, coord, ci = TRUE, ...) {
        forLine <- transform(data, alpha = 1)
        forRibbon <- transform(
            data, 
            alpha = data$alpha, 
            fill = data$colour, 
            colour = NA
        )
        grid::gList(
            if (ci) ggplot2::GeomRibbon$draw_panel(forRibbon, params, coord, ...),
            ggplot2::GeomLine$draw_panel(forLine, params, coord, ...)
        )
    }, 
    draw_key = ggplot2::draw_key_smooth
)

GeomCoverage <- ggplot2::ggproto("GeomCoverage", ggplot2::Geom,
    setup_params = function(data, params) {
        params$type <- params$type
        params$alpha <- params$alpha
        params
    },
    extra_params = c("na.rm", "alpha"),
    required_aes = c("x", "y"), 
    default_aes = ggplot2::aes(
        colour = "black", 
        fill = "grey",
        linewidth = 0.4, 
        linetype = 1, 
        alpha = 1
    ),
    
    draw_group = function(data, params, coord, type, ...) {

        if (!is.null(params$alpha)) {
            data$alpha <- params$alpha
        }
        forArea <- transform(data, ymax = y, ymin = 0, colour = NA)

        grid::gList(
            if (type == 'line') ggplot2::GeomLine$draw_panel(data, params, coord, ...), 
            if (type == 'area') ggplot2::GeomArea$draw_panel(forArea, params, coord, ...)
        )

    }, 
    
    draw_key = function(data, params, type, ...) { 
        if (params$type == 'line') { 
            ggplot2::draw_key_path(data, params)
        } 
        else {
            data <- transform(data, colour = NA)
            ggplot2::draw_key_rect(data, params)
        }
    }
)

#' @rdname ggplot-tidyCoverage
#' @export

geom_aggrcoverage <- function(
    mapping = NULL, 
    data = NULL, 
    ..., 
    unit = c('kb', 'Mb', 'b'), 
    ci = TRUE, 
    grid = FALSE, 
    na.rm = FALSE, 
    show.legend = NA, 
    inherit.aes = TRUE
) {
    m <- ggplot2::aes(x = coord, y = mean, ymin = ci_low, ymax = ci_high, group = interaction(.sample, .feature))
    if (!is.null(mapping)) m <- utils::modifyList(m, mapping)

    unit = match.arg(unit, c('kb', 'Mb', 'b'))
    
    list(
        ggplot2::layer(
            data = data, 
            mapping = m,  
            stat = "identity", 
            geom = GeomAggrCoverage, 
            position = "identity", 
            show.legend = show.legend, 
            inherit.aes = inherit.aes,
            params = list(na.rm = na.rm, ci = ci, ...)
        ), 
        theme_coverage2(grid = grid),
        scale_x_genome(unit = unit)
    )
}

#' @rdname ggplot-tidyCoverage
#' @export

geom_coverage <- function(
    mapping = NULL, 
    data = NULL,
    ..., 
    type = c('area', 'line'), 
    unit = c('kb', 'Mb', 'b'), 
    grid = FALSE, 
    alpha = 0.6,
    na.rm = FALSE, 
    show.legend = NA, 
    inherit.aes = TRUE, 
    raster = TRUE
) {
    m <- ggplot2::aes(x = coord, y = coverage, group = interaction(track, features), fill = track)
    m_line <- ggplot2::aes(x = coord, y = coverage, group = interaction(track, features), color = track)
    if (!is.null(mapping)) m <- utils::modifyList(m, mapping)
    
    unit = match.arg(unit, c('kb', 'Mb', 'b'))
    type <- match.arg(type, c('area', 'line'))

    l <- list(
        ggplot2::layer(
            data = data, 
            mapping = m,  
            stat = "identity", 
            geom = GeomCoverage, 
            position = "identity", 
            show.legend = show.legend, 
            inherit.aes = inherit.aes,
            params = list(na.rm = na.rm, alpha = alpha, type = type, ...)
        ), 
        ggplot2::geom_line(
            data = data,
            mapping = m_line,
            na.rm = na.rm, 
            ...
        ),
        scale_x_genome(unit = unit),
        scale_y_coverage(), 
        theme_coverage(grid = grid), 
        ggplot2::guides(y = ggplot2::guide_axis(cap = "both"))
    ) 
    
    if (raster) {
        l <- l |> ggrastr::rasterize(dpi = 300)
    }
    
    l
}

#' @rdname ggplot-tidyCoverage
#' @export

scale_y_coverage <- function() {
    ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.05)), 
        n.breaks = 3
    )
}

#' @rdname ggplot-tidyCoverage
#' @export

scale_x_genome <- function(unit = c('kb', 'Mb', 'b')) {
    unit = match.arg(unit, c('kb', 'Mb', 'b'))
    scale = dplyr::case_when(
        unit == 'b' ~ 1, 
        unit == 'kb' ~ 1e-3, 
        unit == 'Mb' ~ 1e-6
    )
    ggplot2::scale_x_continuous(
        expand = c(0, 0), 
        labels = scales::unit_format(
            unit = unit, scale = scale, 
            sep = "", 
            big.mark = ""
        )
    )
}

.theme_coverage <- function(
    grid = TRUE, 
    base_size = 11, 
    base_family = "", 
    base_line_size = base_size/22, 
    base_rect_size = base_size/22
) {
    th <- ggplot2::theme_bw(
        base_size = base_size, 
        base_family = base_family, 
        base_line_size = base_line_size, 
        base_rect_size = base_rect_size
    ) 
    if (!grid) th <- th %+replace% ggplot2::theme(
        panel.grid = ggplot2::element_blank(), 
        panel.grid.major = ggplot2::element_blank(), 
        panel.grid.minor = ggplot2::element_blank()
    )
    th <- th %+replace% 
        ggplot2::theme(
            legend.position = 'top', 
            legend.background = ggplot2::element_blank(), 
            legend.key = ggplot2::element_blank(), 
            panel.spacing = unit(8, "pt"), 
            panel.background = ggplot2::element_blank(), 
            strip.background = ggplot2::element_blank(), 
            plot.background = ggplot2::element_blank(), 
            complete = TRUE
        )
    th
}

theme_coverage <- function(
        grid = TRUE, 
        base_size = 11, 
        base_family = "", 
        base_line_size = base_size/22, 
        base_rect_size = base_size/22
) {
    th <- .theme_coverage(
        grid = grid, 
        base_size = base_size, 
        base_family = base_family, 
        base_line_size = base_line_size, 
        base_rect_size = base_rect_size
    ) %+replace% 
        ggplot2::theme(
            #panel.border = ggplot2::element_blank(), 
            axis.line = element_line(color = 'black'), 
            complete = TRUE
        )
    th
}

theme_coverage2 <- function(
    grid = TRUE, 
    base_size = 11, 
    base_family = "", 
    base_line_size = base_size/22, 
    base_rect_size = base_size/22
) {
    th <- .theme_coverage(
        grid = grid, 
        base_size = base_size, 
        base_family = base_family, 
        base_line_size = base_line_size, 
        base_rect_size = base_rect_size
    ) %+replace% 
        ggplot2::theme(
            axis.ticks = ggplot2::element_blank(), 
            complete = TRUE
        )
    th
}