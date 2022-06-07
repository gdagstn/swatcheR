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
