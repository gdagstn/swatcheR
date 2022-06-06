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
#' @details Internal use only. Slightly modified for vectorization.
#'     Originally present in the \code{colorscience} package by Jose Gama.
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
#' @param file_path a character indicating a path to a .JPG or .PNG file. Either "file"
#'     or "link" must be specified.
#' @param link a character indicating a URL to a .JPG or .PNG file. Either "link"
#'     or "file" must be specified.
#' @param m a raster object created by either \code{readJPEG} or \code{readPNG}
#' @param reference_space a reference space array created by a call to \code{makeReferenceSpace()}.
#'    Default is NULL, which internally loads a reference space comprising 4096
#'    colors in DIN99 space.
#' @param keep_extremes logical, should the extreme points of the space
#'    (corresponding to any combination of full R, G, or B channel) be kept?
#'    Default is TRUE.
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
#' @importFrom farver decode_native encode_native
#' @importFrom reshape2 melt
#' @importFrom tools file_ext
#'
#' @export

analyzePictureCol = function(file_path = NULL, link = NULL, m = NULL, reference_space = NULL, keep_extremes = TRUE) {

  if(is.null(file_path) & is.null(link) & is.null(m)) stop("A URL (link), file path (file) or raster (m) must be specified.")
  if(is.null(reference_space)) reference_space = din_4096

  if(is.null(m)){
    if(!is.null(link) & is.null(file_path)) {
      temp = tempfile()
      download.file(link, destfile = temp)
      ext = tolower(tools::file_ext(link))
      if(ext == "jpg" | ext == "jpeg") {
        m = readJPEG(temp, native = TRUE)
      } else if(ext == "png") {
        m = readPNG(temp, native = TRUE)
      }

    } else if(is.null(link) & !is.null(file_path)) {
      ext = tolower(tools::file_ext(file_path))
      if(ext == "jpg" | ext == "jpeg") {
        m = readJPEG(file_path, native = TRUE)
      } else if(ext == "png") {
        m = readPNG(file_path, native = TRUE)
      }
    }
  }

  m2 = melt(matrix(as.numeric(m), nrow = attr(m, "dim")[1], byrow = TRUE))

  m2$hex = decode_native(m2$value)
  m2$value <- NULL
  rgb_values = round(coords(colorspace::hex2RGB(m2$hex)) * 255) + 1

  m2$cluster = reference_space[cbind(rgb_values[,1], rgb_values[,2], rgb_values[,3])]

  if(keep_extremes) {
    for(i in c("#000000", "#FFFFFF", "#FF0000", "#00FF00", "#0000FF","#00FFFF",
               "#FF00FF", "#FFFF00")) {
      m2$cluster[m2$hex == i] = encode_native(i)
    }
  }

  res = table(m2$cluster)
  names(res) = decode_native(names(res))
  return(res)
}

#' Get a summary palette
#'
#' Automatically finds a color swatch from a color analysis
#'
#' @param file_path a character indicating a path to a .JPG or .PNG file.
#'     Either \code{file_path},  \code{link} or  \code{m} must be specified.
#'     Default is NULL.
#' @param link a character indicating a URL to a .JPG or .PNG file. Either
#'     \code{file_path},  \code{link} or  \code{m} must be specified.
#'     Default is NULL.
#' @param m a raster object created by either \code{readJPEG} or \code{readPNG}.
#'     Either \code{file_path},  \code{link} or  \code{m} must be specified.
#'     Default is NULL.
#' @param analysis a named numeric with the results of a color analysis, i.e. the output
#'     of \code{analyzePaintingCol}. Default is NULL.
#' @param reference_space a reference space array created by a call to \code{makeReferenceSpace()}.
#'    Default is NULL, which internally loads a reference space comprising 4096
#'    colors in DIN99 space.
#' @param method a character specifying the method for clustering:
#'     "HC" for hierarchical clustering, "kmeans" for k-means clustering, and
#'     "kmeans_classic" for a legacy implementation. Default is "HC".
#' @param n a numeric, the number of major colors (k-means clusters)
#' @param sub a numeric, the number o minor colors (top abundant colors per cluster)
#' @param ntop a numeric, the number of top colors to use for the initial clustering.
#'     Default is 2000.
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
#' @param keep_extremes logical, passed to \code{analyzePictureCol}
#' @param optimize logical, should the palette be optimized for a minimum
#'     color difference (specified by \code{DE2000_target})? Default is FALSE.
#' @param DE2000_target the minimum DeltaE2000 difference for optimization.
#'
#' @return a character vector of hexadecimal color values of length \code{n} * \code{sub}.
#'
#' @details This function performs several operations on a color analysis result
#'    to isolate a pre-determined number of the most representative colors in the
#'    analysis. First, the analysis results are sorted and the \code{ntop} top colors
#'    are retained; then they are optionally transformed to polarLUV coordinates
#'    (i.e. HCL) and filtered based on the distribution of luminance and/or chroma.
#'    Then, they are transformed to DIN99 space,and clustered using either
#'    hierarchical or k-means clustering, using the user-defined parameter \code{n}
#'    as number of clusters. Another method that does not use the reference space
#'    is "kmeans_classic", which applies k-means directly to the image (after
#'    the RGB values have been transformed to DIN99 space).
#'    Within each cluster, the top \code{sub} colors (as defined by the analysis)
#'    are kept. Optionally, the palette can be optimized by running the same procedure
#'    iteratively until a minimum DeltaE2000 color difference is reached. Optimization
#'    is stopped after 100 unsuccessful attempts, at which point it is advisable to
#'    lower the minimum difference. It is important to note that the bigger the
#'    palette, the lowest the minimum difference between colors will be.
#'    The palette is then optionally ordered according to one of hue, chroma, or
#'    luminance.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace hex2RGB coords
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats kmeans
#' @importFrom methods as
#'
#' @export

getPalette = function(file_path = NULL, link = NULL, m = NULL, analysis = NULL,
                      reference_space = NULL, method = "HC", n = 10, sub = 1,
                      ntop = 2000, filter_luminance = "both", filter_chroma = NULL,
                      filter_sd = 1.5, order = "L", keep_extremes = TRUE,
                      optimize = FALSE, DE2000_target = 4) {


  if(!method %in% c("HC", "kmeans", "kmeans_classic")) stop("Unknown method.")

 if(is.null(analysis)) {
  if(is.null(m)){
    if(is.null(file_path) & is.null(link)) stop("A URL (link), file path (file_path), raster (m) or analysis must be specified.")

    if(!is.null(link) & is.null(file_path)) {
      temp = tempfile()
      download.file(link, destfile = temp)
      ext = tolower(tools::file_ext(link))
      if(ext == "jpg" | ext == "jpeg") {
        m = readJPEG(temp, native = TRUE)
      } else if(ext == "png") {
        m = readPNG(temp, native = TRUE)
      }

    } else if(is.null(link) & !is.null(file_path)) {
      ext = tolower(tools::file_ext(file_path))
      if(ext == "jpg" | ext == "jpeg") {
        m = readJPEG(file_path, native = TRUE)
      } else if(ext == "png") {
        m = readPNG(file_path, native = TRUE)
      }
    }
   }
  }

  if(method == "kmeans_classic") {

    m2 = melt(matrix(as.numeric(m), nrow = attr(m, "dim")[1], byrow = TRUE))
    m2$hex <- hexcols <- decode_native(m2$value)
    labspace = coords(as(hex2RGB(m2$hex), "LAB"))
    din = CIELabtoDIN99mod(labspace[,1], labspace[,2], labspace[,3])

    kk = suppressWarnings(kmeans(din, centers = n))

    m2$cluster = kk$cluster

    if(!is.null(filter_luminance)) {
     luminance_pick = filterColors(hexcols, on = "L",
                                  type = filter_luminance, filter_sd = filter_sd)
    } else luminance_pick = seq_len(length(hexcols))

    if(!is.null(filter_chroma)) {
      chroma_pick = filterColors(hexcols, on = "C",
                                    type = filter_chroma, filter_sd = filter_sd)

    } else chroma_pick = seq_len(length(hexcols))

    hexcols = hexcols[intersect(luminance_pick, chroma_pick)]
    m2 = m2[m2$hex %in% hexcols,]
    m2$pixels = table(m2$hex)[m2$hex]
    m2 = m2[order(m2$pixels, decreasing = TRUE),]
    m2_list = split(m2, m2$cluster)

    pal_kdin = unlist(lapply(m2_list, function(x) x[seq_len(sub),"hex"]), use.names = FALSE)

  } else if(method != "kmeans_classic") {

    if(is.null(analysis)) analysis = analyzePictureCol(m = m, reference_space = reference_space, keep_extremes = keep_extremes)

    sorted = sort(analysis, decreasing = TRUE)[seq_len(min(ntop, sum(analysis > 0)))]

    hexcols = names(sorted)

    if(!is.null(filter_luminance)) {
      luminance_pick = filterColors(hexcols, on = "L",
                                    type = filter_luminance, filter_sd = filter_sd)
    } else luminance_pick = seq_len(length(hexcols))

    if(!is.null(filter_chroma)) {
      chroma_pick = filterColors(hexcols, on = "C",
                                 type = filter_chroma, filter_sd = filter_sd)

    } else chroma_pick = seq_len(length(hexcols))

    hexcols = hexcols[intersect(luminance_pick, chroma_pick)]

    sorted = sorted[hexcols]
    sorted = sort(sorted, decreasing = TRUE)

    pal_kdin = clusterKDIN(sorted = sorted, n = n, sub = sub, method = method)
  }

    if(optimize) {
      counter = 0
      pal_kdin_temp = pal_kdin
      to_optimize = min(getPaletteDistances(pal_kdin_temp)$DE2000)
      while(to_optimize < DE2000_target) {
        pal_kdin_temp = clusterKDIN(sorted = sorted, n = n, sub = sub, method = method)
        to_optimize = min(getPaletteDistances(pal_kdin_temp)$DE2000)
        counter = counter + 1
        if(counter > 99) {
          message(paste0("Failed to optimize in ", counter, " iterations. Try reducing the DE2000 target."))
          break()
         }
        }
        if(counter <= 99) message(paste0("Optimized in ", counter, " iterations."))
        pal_kdin = pal_kdin_temp
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
#' @param file_path a character indicating a path to a .JPG or .PNG file. Either "file"
#'     or "link" must be specified.
#' @param link a character indicating a URL to a .JPG or .PNG file. Either "link"
#'     or "file" must be specified.
#' @param m a raster object created by either \code{readJPEG} or \code{readPNG}
#' @param reference_space a reference space array created by a call to \code{makeReferenceSpace()}
#' @param method a character specifying the method for clustering:
#'     "HC" for hierarchical clustering, "kmeans" for k-means clustering, and
#'     "kmeans_classic" for a legacy implementation. Default is "HC".
#' @param n a numeric, the number of major colors (clusters)
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
#' @param optimize logical, should the palette be optimized for maximal difference?
#'     Default is FALSE
#' @param DE2000_target numeric, the target minimum color difference for optimization.
#'     Default is 4. Only relevant if \code{optimize} is TRUE.
#' @param bg_color the background color (as named color, hexadecimal or rgb).
#'     Default is "white".
#' @param keep_extremes logical, passed to \code{analyzePictureCol} if needed.
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
#' @importFrom ggplot2 ggplot geom_rect geom_text theme_void annotation_raster theme_void theme ggtitle aes_string element_rect coord_fixed xlim ylim element_text
#'
#' @export

plotWithPal = function(file_path = NULL, link = NULL, m = NULL, reference_space = NULL,
                       method = "HC", n = 20, sub = 1, filter_luminance = "both",
                       filter_chroma = NULL, filter_sd = 1.5, order = "L",
                       optimize = FALSE, DE2000_target = 4, bg_color = "white",
                       keep_extremes = TRUE, title = NULL){

  if(is.null(m)){

    if(is.null(file_path) & is.null(link)) stop("A URL (link), file path (file_path) or raster (m) must be specified.")

  if(!is.null(link)) {
    temp = tempfile()
    download.file(link, destfile = temp)
    path = temp
    ext = tolower(file_ext(link))
  } else if(!is.null(file_path)) {
    path = file_path
    ext = tolower(file_ext(path))
  }

  if(ext == "jpg" | ext == "jpeg") {
    m = readJPEG(path, native = TRUE)
  } else if(ext == "png") {
    m = readPNG(path, native = TRUE)
  }
 }
  pal = getPalette(m = m, reference_space = reference_space,
                   method = method, n = n, sub = sub,
                   filter_luminance = filter_luminance,
                   filter_chroma = filter_chroma, filter_sd = filter_sd,
                   order = order, optimize = optimize,
                   DE2000_target = DE2000_target)

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

  if(!is.null(title)) p = p +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)

}

#' Get palette distances
#'
#' Calculate pairwise differences between colors in a palette
#'
#' @param palette a vector of colors as hexadecimal values
#'
#' @return a \code{data.frame} with color pairs and their DeltaE2000 color difference
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace hex2RGB
#' @importFrom methods as

getPaletteDistances = function(palette) {

  paletteLAB = as(hex2RGB(palette), "LAB")

  dists = expand.grid(seq_len(length(palette)), seq_len(length(palette)))
  dists = dists[dists[,1] != dists[,2],]
  dists$Col1 = palette[dists[,1]]
  dists$Col2 = palette[dists[,2]]
  dists$DE2000 = deltaE2000mod(coords(paletteLAB)[dists$Var1,], coords(paletteLAB)[dists$Var2,])

  return(dists)
}

#' Palette hierarchical clustering
#'
#' Applies hierarchical clustering to a palette returning clusters
#'
#' @param palette a vector of colors as hexadecimal values
#' @param n a numeric, desired number of clusters
#'
#' @return a \code{data.frame} with colors and clusters.
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom reshape2 acast
#' @importFrom stats hclust dist cutree as.dist

hClusterPalette = function(palette, n) {

  df = getPaletteDistances(palette)
  distmat = as.dist(acast(df, formula = Var1 ~ Var2, value.var = "DE2000"))
  hc = hclust(distmat)
  clusters = cutree(hc, k = n)

  clusters_df = data.frame("col" = palette[as.numeric(names(clusters))], "cluster" = clusters)

  return(clusters_df)
}

#' DeltaE200 color difference
#'
#' Calculates the Delta E (CIE 2000) color difference between two matrices in CIELab space
#'
#' @param Lab1 matrix with L, a, b coordinates for one set of colors
#' @param Lab2 matrix with L, a, b coordinates for another set of colors
#'
#' @return a vector of DeltaE2000 color distances
#'
#' @details Internal use only. Slightly modified for vectorisation.
#'     Originally present in the \code{colorscience} package by Jose Gama.
#'
#' @author Jose Gama, modified by Giuseppe D'Agostino
#'
#' @references Bruce Justin Lindbloom, 2013 Color Calculator http://www.brucelindbloom.com

deltaE2000mod <- function (Lab1, Lab2) {
  kL <- 1
  kC <- 1
  kH <- 1
  lBarPrime <- 0.5 * (Lab1[,1] + Lab2[,1])
  c1 <- sqrt(Lab1[,2] * Lab1[,2] + Lab1[,3] * Lab1[,3])
  c2 <- sqrt(Lab2[,2] * Lab2[,2] + Lab2[,3] * Lab2[,3])
  cBar <- 0.5 * (c1 + c2)
  cBar7 <- cBar^7
  g <- 0.5 * (1 - sqrt(cBar7/(cBar7 + 6103515625)))
  a1Prime <- Lab1[,2] * (1 + g)
  a2Prime <- Lab2[,2] * (1 + g)
  c1Prime <- sqrt(a1Prime * a1Prime + Lab1[,3] * Lab1[,3])
  c2Prime <- sqrt(a2Prime * a2Prime + Lab2[,3] * Lab2[,3])
  cBarPrime <- 0.5 * (c1Prime + c2Prime)
  h1Prime <- (atan2(Lab1[,3], a1Prime) * 180)/pi
  h1Prime[h1Prime < 0] <- h1Prime[h1Prime < 0] + 360
  h2Prime <- (atan2(Lab2[,3], a2Prime) * 180)/pi
  h2Prime[h2Prime < 0] <- h2Prime[h2Prime < 0] + 360

  hBarPrime <- abs(h1Prime - h2Prime)
  hBarPrime[hBarPrime > 180] <-  0.5 * (h1Prime[hBarPrime > 180] + h2Prime[hBarPrime > 180] + 360)
  hBarPrime[hBarPrime <= 180] <- 0.5 * (h1Prime[hBarPrime <= 180] + h2Prime[hBarPrime <= 180])

  t <- 1 - 0.17 * cos(pi * (hBarPrime - 30)/180) + 0.24 * cos(pi * (2 * hBarPrime)/180) +
    0.32 * cos(pi * (3 * hBarPrime + 6)/180) - 0.2 * cos(pi * (4 * hBarPrime - 63)/180)

  dhPrime = abs(h2Prime - h1Prime)
  dhPrime[dhPrime <= 180] <- h2Prime[dhPrime <= 180] - h1Prime[dhPrime <= 180]
  dhPrime[dhPrime > 180 & h2Prime <= h1Prime] <- h2Prime[dhPrime > 180 & h2Prime <= h1Prime]  -
    h1Prime[dhPrime > 180 & h2Prime <= h1Prime]  + 360
  dhPrime[dhPrime > 180 & h2Prime > h1Prime] <- h2Prime[dhPrime > 180 & h2Prime > h1Prime] -
    h1Prime[dhPrime > 180 & h2Prime > h1Prime] - 360

  dLPrime <- Lab2[,1] - Lab1[,1]
  dCPrime <- c2Prime - c1Prime
  dHPrime <- 2 * sqrt(c1Prime * c2Prime) * sin(pi * (0.5 *
                                                       dhPrime)/180)
  sL <- 1 + ((0.015 * (lBarPrime - 50) * (lBarPrime - 50))/sqrt(20 +
                                                                  (lBarPrime - 50) * (lBarPrime - 50)))

  sC <- 1 + 0.045 * cBarPrime
  sH <- 1 + 0.015 * cBarPrime * t
  dTheta <- 30 * exp(-((hBarPrime - 275)/25) * ((hBarPrime -
                                                   275)/25))
  cBarPrime7 <- cBarPrime * cBarPrime * cBarPrime * cBarPrime *
    cBarPrime * cBarPrime * cBarPrime
  rC <- sqrt(cBarPrime7/(cBarPrime7 + 6103515625))
  rT <- -2 * rC * sin(pi * (2 * dTheta)/180)
  res = sqrt((dLPrime/(kL * sL)) * (dLPrime/(kL * sL)) +
               (dCPrime/(kC * sC)) * (dCPrime/(kC * sC)) +
               (dHPrime/(kH * sH)) * (dHPrime/(kH * sH)) +
               (dCPrime/(kC * sC)) * (dHPrime/(kH * sH)) * rT)
  return(res)
}


#' Get kNN graph
#'
#' Builds a k-nearest neighbor graph from a set of points
#'
#' @param coords matrix with point coordinates
#' @param k number of nearest neighbors
#'
#' @return an \code{igraph} graph object with the k-NN graph
#'
#' @details Internal use only. Nothing fancy, just a brute force
#'     barebones implementation of k-NN.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom igraph simplify graph_from_data_frame
#' @importFrom stats dist


getKNN <- function(coords, k = 2) {

  distmat = as.matrix(dist(coords))
  diag(distmat) = NA

  index = t(apply(distmat, 1, function(x) order(x)[seq_len(k)]))
  indices = list()
  for(i in seq_len(nrow(index))) {
    indices[[i]] <- list()
    for(j in seq_len(ncol(index))){
      indices[[i]][[j]] <- c(index[i,1], index[i,j])
    }
  }

  index_nn = as.data.frame(do.call(rbind, lapply(indices, function(x) do.call(rbind, x))))
  g = simplify(graph_from_data_frame(index_nn, directed = FALSE))
  return(g)
}

#' Get continuous palettes
#'
#' Generates continuous palettes from a categorical palette
#'
#' @param pal a hex code vector with the palette
#' @param n a numeric, the final number of palette colors. Default is 10
#' @param tries numeric, how many times should the algorithm run?
#' @param plot logical, should the resulting palettes be plotted?
#'
#' @return a nested list of length \code{tries}, each element containing:
#' \itemize{
#'    \item{the original palette in DIN99 space (\code{original_space})}
#'    \item{the palette space after clustering (\code{palette_space})}
#'    \item{the principal curve (\code{principal_curve})}
#'    \item{the curve colors (\code{principal_curve_rgb})}
#'    \item{the final palette (\code{palette})}
#'    }
#'
#' @details todo
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom igraph subgraph components shortest_paths mst distances
#' @importFrom stats kmeans
#' @importFrom princurve principal_curve
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom methods as
#' @importFrom colorspace LAB coords swatchplot

getContinuousPalette <- function(pal, n = 10, tries = 1, plot = TRUE) {

  palettes = list()

  for(i in seq_len(tries)){

    cielab = coords(as(hex2RGB(pal), "LAB"))
    din = CIELabtoDIN99mod(cielab[,1], cielab[,2], cielab[,3])
    dink = kmeans(din, centers = round(sqrt(length(pal))))

    g = getKNN(dink$centers, k = 2)
    comp = which.max(table(components(g)$membership))
    vkeep = names(components(g)$membership[components(g)$membership == comp])
    g = subgraph(g, vids = vkeep)

    startfrom = which.max(rowSums(distances(g)))
    path_chosen = which.max(lengths(shortest_paths(mst(g), from = startfrom, algorithm = "bellman-ford")$vpath))
    path_vs = as.numeric(names(shortest_paths(mst(g), from = startfrom, algorithm = "bellman-ford")$vpath[[path_chosen]]))

    keep = din[dink$cluster %in% as.numeric(path_vs),]
    pcur = principal_curve(keep)
    pcur_colors = coords(as(LAB(DIN99toCIELabmod(pcur$s[,1], pcur$s[,2], pcur$s[,3])), "sRGB"))
    pcur_colors[pcur_colors > 1] = 1
    pcur_colors[pcur_colors < 0] = 0
    pcur_rgb = rgb(pcur_colors)[pcur$ord]

    palettes[[i]] = list("original_space" = din,  "palette_space" = keep, "principal_curve" = pcur, "principal_curve_rgb" = pcur_rgb,
                         "palette" = colorRampPalette(pcur_rgb)(n))
  }

  if(plot) swatchplot(lapply(palettes, function(x) unlist(x$palette)))
  return(palettes)
}


#' K clusters in DIN99 space
#'
#' Clusters a color palette in DIN99 color space returning k clusters
#'
#' @param sorted picture analysis from \code{analyzePictureCol}
#' @param n a numeric, the final number of clusters
#' @param sub a numeric, most abundant colors per cluster
#' @param method a character, one of "kmeans" or "HC"
#'
#' @return a palette with most abundant and separated colors
#'
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace coords hex2RGB
#' @importFrom stats kmeans
#' @importFrom grDevices col2rgb
#' @importFrom methods as

clusterKDIN <- function(sorted, n, sub, method){

  hexcols = names(sorted)

  if(method == "kmeans") {
    cielab = coords(as(hex2RGB(hexcols), "LAB"))

    din =  CIELabtoDIN99mod(L = cielab[,1], a = cielab[,2], b = cielab[,3])
    kdin = kmeans(din, centers = min(c(n, nrow(din) - 1)), nstart = 5)
    ks_sorted = data.frame("pixels" = as.numeric(sorted),
                           "cluster" = kdin$cluster, row.names = names(sorted))

  } else if(method == "HC"){
    ks_sorted = hClusterPalette(palette = hexcols, n = n)
    ks_sorted$pixels =  sorted[ks_sorted$col]
    rownames(ks_sorted) = ks_sorted$col
    ks_sorted$col = NULL
    ks_sorted = ks_sorted[,2:1]
  }

  pal_kdin = as.character(unlist(lapply(split(ks_sorted, ks_sorted$cluster), function(x) {
    if(nrow(x) >= sub) {
      return(rownames(x)[seq_len(sub)])
    } else {
      return(rownames(x))}}), use.names = FALSE))

  if(length(pal_kdin) < n * sub) {
    pal_kdin = c(pal_kdin, rep("#00000000", (n * sub)-length(pal_kdin)))
  }

  return(pal_kdin)
}


#' Filter colors
#'
#' Filters colors based on Luminace or Chroma
#'
#' @param colors hex code vector of colors
#' @param on character, one of "L" (luminance) or "C" (chroma)
#' @param type character, one of "dark" or "bright" for "L",
#'     "high" or "low" for "C", and "both" for either value of \code{on}
#' @param filter_sd numeric, the number of SD/MAD for filtering
#'
#' @return a vector of indices for filtering colors
#'
#' @details internal use only.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom colorspace hex2RGB
#' @importFrom methods as
#' @importFrom stats sd mad median

filterColors <- function(colors, on, type, filter_sd) {

  pal_RGB = hex2RGB(colors)

  if(on == "L"){
    luminance = coords(as(pal_RGB, "polarLUV"))[,on]
    if(type == "both") {
      pick = which(luminance <= (mean(luminance) + (filter_sd * sd(luminance))) & luminance >= (mean(luminance) - (filter_sd * sd(luminance))))
    } else if(type == "bright") {
      pick = which(luminance >= (mean(luminance) - (filter_sd * sd(luminance))))
    } else if(type == "dark") {
      pick = which(luminance <= (mean(luminance) + (filter_sd * sd(luminance))))
    }
  }

  if(on == "C"){
    chroma = coords(as(pal_RGB, "polarLUV"))[,on]

    if(type == "both") {
      pick = which(chroma <= (median(chroma) + (filter_sd * mad(chroma))) & chroma >= (median(chroma) - (filter_sd * mad(chroma))))
    } else if(type == "high") {
      pick = which(chroma >= (median(chroma) - (filter_sd * mad(chroma))))
    } else if(type == "low") {
      pick = which(chroma <= (median(chroma) + (filter_sd * mad(chroma))))
    }
  }

  return(pick)
}
