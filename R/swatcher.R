#' Make reference space
#'
#' Build a summarized color reference space
#'
#' @param k a numeric indicating the final number of summary colors. Default is 1000
#' @param space a character, one of "LAB" or "DIN99"
#' @param verbose logical, should the function output information about the progress? Default is FALSE
#'
#' @return a cube 3D \code{array} of side 256 in which every value corresponds to a
#'     summarized color for that RGB coordinate. Values are encoded as native integers.
#'
#' @details This function creates the reference space to assign single RGB colors
#'     (defined by a triplet of R, G, B coordinates between 0 and 255) to one of
#'     a number (k) of clusters. Clusters are defined via k-means within one of
#'     two perceptual color spaces, CIELab or DIN99. The function builds the RGB
#'     cube first, then transforms it to either CIELab or DIN99 coordinates, performs
#'     k-means clustering and takes cluster centroid coordinates as summary colors.
#'     Summary colors are then assigned back to the original RGB cube. For a particular
#'     RGB triplet the corresponding position in the array is matched with the summary
#'     color (centroid of the clustering it was included in).
#'
#' @note Since the number \code{k} of clusters is usually high (1000-3000) and
#'     distances are regularly spaced, the Hartigan-Wong algorithm may not converge
#'     in Quick-Transfer stage. This will throw a warning, but it does not seem to
#'     affect the results for this application.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace RGB LAB coords
#' @importFrom farver encode_native
#' @importFrom reshape2 acast
#' @importFrom stats kmeans
#' @importFrom methods as
#'
#' @export

makeReferenceSpace = function(k = 1000, space = "LAB", verbose = FALSE){

  if(k <= 1) stop("k must be greater than 1")
  if(!space %in% c("LAB", "DIN99")) stop("space must be one of \"LAB\" or \"DIN99\"")

  rgbcube = expand.grid(0:255, 0:255, 0:255)
  rgbcube = RGB(rgbcube[,1]/255, rgbcube[,2]/255, rgbcube[,3]/255)

  if(verbose) message(paste0("K-means clustering ", space, " space with k = ", k, "..."))

  if(space == "LAB") {
    space = coords(as(rgbcube, "LAB"))
    kmeans_space = kmeans(space, centers = k)
    kmeans_centers = LAB(kmeans_space$centers[,1], kmeans_space$centers[,2], kmeans_space$centers[,3])

  } else if(space == "DIN99") {
    lab = coords(as(rgbcube, "LAB"))
    space = CIELabtoDIN99mod(lab[,1], lab[,2], lab[,3])
    kmeans_space = kmeans(space, centers = k)
    kmeans_centers = LAB(DIN99toCIELabmod(kmeans_space$centers[,1], kmeans_space$centers[,2], kmeans_space$centers[,3]))
  }

  kmeans_centers_rgb = coords(as(kmeans_centers, "RGB"))

  if(verbose) message(paste0("Creating encoded array..."))
  kmeans_cols = as.numeric(encode_native(rgb(kmeans_centers_rgb[,1], kmeans_centers_rgb[,2], kmeans_centers_rgb[,3])))

  rgb_with_cols = data.frame("R" = coords(rgbcube)[,1],
                             "G" = coords(rgbcube)[,2],
                             "B" = coords(rgbcube)[,3],
                             "col" = kmeans_cols[kmeans_space$cluster])

  final_arr = acast(rgb_with_cols, formula = R ~ G ~ B, value.var = "col")

  if(verbose) message(paste0("Done."))
  return(final_arr)
}

#' CIELab to DIN99 transformation
#'
#' Transform L, a, b, to L99o, a99o, b99o
#'
#' @param L numeric vector of L values
#' @param a numeric vector of a values
#' @param b numeric vector of b values
#'
#' @return a matrix with the same number of rows as the length of the vectors in
#'     DIN99 coordinates.
#'
#' @details Internal use only. Slightly modified to take vectors and not single
#'     values as input. Originally present in the \code{colorscience} package
#'     by Jose Gama.
#'
#' @author Jose Gama, modified by Giuseppe D'Agostino
#'
#' @references CIELAB to DIN99 coordinates, 2014 http://de.wikipedia.org/w/index.php?title=Diskussion:DIN99-Farbraum

CIELabtoDIN99mod <- function (L, a, b) {

  if(mean(length(L), length(a), length(b)) != length(L)) stop("L, a, and b must have the same length")

  kE <- 1
  kCH <- 1
  ang <- 2 * pi/360 * 26
  L99f <- 100/log(139/100)
  L99o <- L99f/kE * log(1 + 0.0039 * L)
  eo <- a * cos(ang) + b * sin(ang)
  fo <- 0.83 * (b * cos(ang) - a * sin(ang))
  Go <- sqrt(eo^2 + fo^2)
  C99o <- log(1 + 0.075 * Go)/(0.0435 * kCH * kE)
  heofo <- atan2(fo, eo)
  h99o <- heofo + ang
  a99o <- C99o * cos(h99o)
  b99o <- C99o * sin(h99o)
  cbind(L99o, a99o, b99o)
}

#' DIN99 to CIELab transformation
#'
#' Transform L99o, a99o, b99o to L, a, b
#'
#' @param L99o numeric vector of L99o values
#' @param a99o numeric vector of a99o values
#' @param b99o numeric vector of b99o values
#'
#' @return a matrix with the same number of rows as the length of the vectors in
#'     LAB coordinates.
#'
#' @details Internal use only. Slightly modified to take vectors and not single
#'     values as input. Originally present in the \code{colorscience} package
#'     by Jose Gama.
#'
#' @author Jose Gama, modified by Giuseppe D'Agostino
#'
#' @references DIN99 to CIELAB coordinates, 2014 http://de.wikipedia.org/w/index.php?title=Diskussion:DIN99-Farbraum

DIN99toCIELabmod <- function (L99o, a99o, b99o) {

  if(mean(length(L99o), length(a99o), length(b99o)) != length(L99o)) stop("L99o, a99o, and b99o must have the same length")

  kE <- 1
  kCH <- 1
  ang <- 2 * pi/360 * 26
  L99f <- 100/log(139/100)
  L <- (exp(L99o * kE/L99f) - 1)/0.0039
  h99ef <- atan2(b99o, a99o)
  heofo <- h99ef - ang
  C99 <- sqrt(a99o^2 + b99o^2)
  G <- (exp(0.0435 * kE * kCH * C99) - 1)/0.075
  e <- G * cos(heofo)
  f <- G * sin(heofo)
  a <- (e * cos(ang) - (f/0.83) * sin(ang))
  b <- (e * sin(ang) + (f/0.83) * cos(ang))
  cbind(L = L, a = a, b = b)
}

#' Analyze picture colors
#'
#' Quantifies summary colors in a picture
#'
#' @param file a character indicating a path to a .JPG or .PNG file. Either "file"
#'     or "link" must be specified.
#' @param link a character indicating a URL to a .JPG or .PNG file. Either "link"
#'     or "file" must be specified.
#' @param reference_space a reference space array created by a call to \code{makeReferenceSpace()}.
#'    Default is NULL, which internally loads a reference space comprising 2000
#'    colors in DIN99 space.
#'
#' @return a named numeric vector of summary color quantification (number of
#'     pixels with that summary color), whose names are the hexadecimal values of
#'     colors.
#'
#' @details This function decomposes a .JPG or .PNG file into component RGB values,
#'     assigns each of these values to a summary color as indicated in the reference
#'     space, and quantifies the number of pixels for each summary color. The result
#'     is unsorted.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom jpeg readJPEG
#' @importFrom png readPNG
#' @importFrom utils download.file
#' @importFrom farver decode_native
#' @importFrom reshape2 melt
#' @importFrom tools file_ext
#'
#' @export


analyzePictureCol = function(file = NULL, link = NULL, reference_space = NULL) {

  if(is.null(file) & is.null(link)) stop("A URL (link) or file path (file) must be specified.")
  if(is.null(reference_space)) reference_space = din_2000

  if(!is.null(link)) {
    temp = tempfile()
    download.file(link, destfile = temp)
    ext = tolower(file_ext(temp))
    if(ext == "jpg" | ext == "jpeg") {
      m = readJPEG(temp, native = TRUE)
    } else if(ext == "png") {
      m = readPNG(temp, native = TRUE)
    }

  } else if(is.null(link) & !is.null(file)) {
    ext = tolower(file_ext(file))
    if(ext == "jpg" | ext == "jpeg") {
      m = readJPEG(file, native = TRUE)
    } else if(ext == "png") {
      m = readPNG(file, native = TRUE)
    }
  }

  m2 = melt(matrix(as.numeric(m), nrow = attr(m, "dim")[1], byrow = TRUE))

  m2$hex = decode_native(m2$value)
  m2$value <- NULL
  rgb_values = round(coords(colorspace::hex2RGB(m2$hex)) * 255) + 1

  m2$cluster = reference_space[cbind(rgb_values[,1], rgb_values[,2], rgb_values[,3])]
  res = table(m2$cluster)
  names(res) = decode_native(names(res))
  return(res)
}



#' Get a summary palette
#'
#' Automatically finds a color palette from a color analysis
#'
#' @param analysis a named numeric with the results of a color analysis, i.e. the output
#'     of \code{analyzePaintingCol}
#' @param n a numeric, the number of major colors (k-means clusters)
#' @param sub a numeric, the number o minor colors (top abundant colors per cluster)
#' @param filter_luminance character, one of "dark", "bright", or "both". Removes
#'     colors within a certain number of standard deviations from the mean of the
#'     luminance values. "dark" keeps darker colors, "bright" keeps brighter colors,
#'     "both" applies both filters. Default is "both".
#' @param filter_chroma character, one of "low", "high", or "both". Removes
#'     colors within a certain number of median absolute deviations from the median
#'     of the chroma values. "low" keeps less saturated colors, "high" keeps
#'     more saturated colors, "both" applies both filters. Default is NULL.
#' @param filter_sd numeric, the number of standard deviations and/or median absolute
#'     deviations above/below which filters are applied. Default is 1.5.
#' @param order character, the dimension of HCL space along which the final palette
#'     should be ordered. One of "H" (hue), "C" (chroma), "L" (luminance). The
#'     final order has only aesthetic purposes and does not change the palette.
#'     Default is "L".
#'
#' @return a character vector of hexadecimal color values of length \code{n} * \code{sub}.
#'
#' @details This function performs several operations on a color analysis result
#'    to isolate a pre-determined number of the most representative colors in the
#'    analysis. First, the analysis results are sorted and the top 500 colors are retained;
#'    then they are optionally transformed to polarLUV coordinates (i.e. HCL) and
#'    filtered based on the distribution of luminance and/or chroma. Then, they
#'    are transformed to DIN99 space,and clustered using k-means clustering, using
#'    the user-defined parameter \code{n} as number of clusters. Within each
#'    cluster, the top \code{sub} colors (as defined by the analysis) are kept.
#'    The palette is then optionally ordered according to one of hue, chroma, or
#'    luminance.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace hex2RGB LAB RGB polarLUV coords
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats kmeans mad median sd
#' @importFrom methods as
#'
#' @export

getPalette = function(analysis, n = 10, sub = 1, filter_luminance = "both",
                      filter_chroma = NULL, filter_sd = 1.5, order = "L") {

  sorted = sort(analysis, decreasing = TRUE)[seq_len(min(500, sum(analysis > 0)))]

  hexcols = names(sorted)

  if(!is.null(filter_luminance)) {

    pal_RGB = hex2RGB(hexcols)
    luminance = coords(as(pal_RGB, "polarLUV"))[,"L"]

    if(filter_luminance == "both") {
      luminance_pick = which(luminance <= (mean(luminance) + (filter_sd * sd(luminance))) & luminance >= (mean(luminance) - (filter_sd * sd(luminance))))
    } else if(filter_luminance == "bright") {
      luminance_pick = which(luminance >= (mean(luminance) - (filter_sd * sd(luminance))))
    } else if(filter_luminance == "dark") {
      luminance_pick = which(luminance <= (mean(luminance) + (filter_sd * sd(luminance))))
    }
  } else luminance_pick = seq_len(length(hexcols))

  if(!is.null(filter_chroma)) {

    pal_RGB = hex2RGB(hexcols)
    chroma = coords(as(pal_RGB, "polarLUV"))[,"C"]

    if(filter_chroma == "both") {
      chroma_pick = which(chroma <= (median(chroma) + (filter_sd * mad(chroma))) & chroma >= (median(chroma) - (filter_sd * mad(chroma))))
    } else if(filter_chroma == "high") {
      chroma_pick = which(chroma >= (median(chroma) - (filter_sd * mad(chroma))))
    } else if(filter_chroma == "low") {
      chroma_pick = which(chroma <= (median(chroma) + (filter_sd * mad(chroma))))
    }
  } else chroma_pick = seq_len(length(hexcols))

  hexcols = hexcols[intersect(luminance_pick, chroma_pick)]

  sorted = sorted[hexcols]

  sorted = sort(sorted, decreasing = TRUE)

  cielab = coords(as(RGB(t(col2rgb(hexcols))[,1]/255, t(col2rgb(hexcols))[,2]/255, t(col2rgb(hexcols))[,3]/255), "LAB"))

  din =  CIELabtoDIN99mod(L = cielab[,1], a = cielab[,2], b = cielab[,3])

  kdin = kmeans(din, centers = min(c(n, nrow(din) - 1)), nstart = 100)

  ks_sorted = data.frame("pixels" = as.numeric(sorted), "cluster" = kdin$cluster, row.names = names(sorted))

  pal_kdin = as.character(unlist(lapply(split(ks_sorted, ks_sorted$cluster), function(x) {
    if(nrow(x) >= sub) {
      return(rownames(x)[seq_len(sub)])
    } else {
      return(rownames(x))}}), use.names = FALSE))

  if(length(pal_kdin) < n * sub) {
    pal_kdin = c(pal_kdin, rep("#00000000", (n * sub)-length(pal_kdin)))
  }

  if(!is.null(order)) {
    colz = hex2RGB(pal_kdin)
    luv = as(colz, "polarLUV")
    pal_final = rgb(coords(colz[order(coords(luv)[,order], decreasing = TRUE),]))
  } else {
    pal_final = pal_kdin
  }

  return(unique(pal_final))
}


#' Plot a picture with its palette
#'
#' Plots a picture and the resulting summary palette side by side
#'
#' @param file a character indicating a path to a .JPG or .PNG file. Either "file"
#'     or "link" must be specified.
#' @param link a character indicating a URL to a .JPG or .PNG file. Either "link"
#'     or "file" must be specified.
#' @param reference_space a reference space array created by a call to \code{makeReferenceSpace()}
#' @param n a numeric, the number of major colors (k-means clusters)
#' @param sub a numeric, the number o minor colors (top abundant colors per cluster)
#' @param filter_luminance character, one of "dark", "bright", or "both". Removes
#'     colors within a certain number of standard deviations from the mean of the
#'     luminance values. "dark" keeps darker colors, "bright" keeps brighter colors,
#'     "both" applies both filters. Default is "both".
#' @param filter_chroma character, one of "low", "high", or "both". Removes
#'     colors within a certain number of median absolute deviations from the median
#'     of the chroma values. "low" keeps less saturated colors, "high" keeps
#'     more saturated colors, "both" applies both filters. Default is NULL.
#' @param filter_sd numeric, the number of standard deviations and/or median absolute
#'     deviations above/below which filters are applied. Default is 1.5.
#' @param order character, the dimension of HCL space along which the final palette
#'     should be ordered. One of "H" (hue), "C" (chroma), "L" (luminance). The
#'     final order has only aesthetic purposes and does not change the palette.
#'     Default is "L".
#' @param bg_color the background color (as named color, hexadecimal or rgb)
#' @param title character, the title of the plot.
#'
#' @return a \code{ggplot2} plot with the picture and the corresponding palette.
#'
#' @details This function is a wrapper around ]code{analyzePaintingCol} and
#'     \code{getPalette} that allows the user to select a picture, extract the
#'     palette and plot them side by side.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace hex2RGB LAB RGB polarLUV coords
#' @importFrom grDevices col2rgb rgb
#' @importFrom methods as
#' @importFrom jpeg readJPEG
#' @importFrom png readPNG
#' @importFrom tools file_ext
#' @importFrom ggplot2 ggplot geom_rect geom_text theme_void annotation_raster theme_void theme ggtitle aes_string element_rect coord_fixed xlim ylim
#'
#' @export

plotWithPal = function(file = NULL, link = NULL, reference_space = NULL,
                       n = 20, sub = 1, filter_luminance = "both",
                       filter_chroma = NULL, filter_sd = 1.5, order = "L",
                       bg_color = "black", title = NULL){

  if(!is.null(link)) {
    temp = tempfile()
    download.file(link, destfile = temp)
    path = temp
    ext = tolower(file_ext(link))
  } else if(!is.null(file)) {
    path = file
    ext = tolower(file_ext(path))
  }

  pal = getPalette(analysis = analyzePictureCol(path, reference_space = reference_space),
                   n = n, sub = sub, filter_luminance = filter_luminance,
                   filter_chroma = filter_chroma, filter_sd = filter_sd,
                   order = order)

  if(ext == "jpg" | ext == "jpeg") {
    m = readJPEG(path, native = TRUE)
  } else if(ext == "png") {
    m = readPNG(path, native = TRUE)
  }

  res = dim(m)[2:1]

  rect_ys = seq(0, res[2], length.out = length(pal) + 1)
  bg_lum = coords(polarLUV(t(col2rgb(bg_color))))[,"L"]
  if(bg_lum <= 128) textcol = "white" else textcol = "black"

  p = ggplot(data = data.frame("x" = 1, "y" = 1), aes_string(x = "x", y = "y")) +
    xlim(1 -diff(rect_ys)[1]*3 , res[1]+ res[1]/20 + diff(rect_ys)[1]*4) +
    ylim(1, res[2]) +
    annotation_raster(m, xmin = 1, ymin = 1, xmax = res[1], ymax = res[2])

  for(i in seq_len(length(pal))) {
    p = p + geom_rect(xmin = res[1] + res[1]/20,
                      xmax = res[1] + res[1]/20 + diff(rect_ys)[1],
                      ymin = rect_ys[i],
                      ymax = rect_ys[i + 1],
                      fill = pal[i], color = "black")

    p = p + geom_text(x = res[1] + res[1]/20 + diff(rect_ys)[1]*1.2,
                      y = rect_ys[i] + (res[2]/(length(pal)*2)),
                      label = pal[i], size = 2.5, color = textcol,
                      hjust = "left")
  }
  p = p +
    coord_fixed() +
    theme_void() +
    theme(plot.background = element_rect(fill = bg_color, color = bg_color))

  if(!is.null(title)) p = p + ggtitle(title)

  return(p)

}

