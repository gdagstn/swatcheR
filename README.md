<img src="https://user-images.githubusercontent.com/21171362/171850532-ded8fd04-6233-4f68-b278-2a23bf1ee495.png" align="right" alt="" width="150" />

# swatcheR
Fast automatic selection of color palettes from images

# Installation

You can install `swatcher` through `devtools`:

```
devtools::install_github("gdagstn/swatcheR")
```

# Motivation and methodology
This package lets you automatically and quickly select a color palette (or **swatch**) from a picture in .jpg or .png format. It is based on k-means clustering, which means that there will be some variability in the results unless you set a seed for the random number generator. 

The idea is not new: the `paletter` [package](https://github.com/AndreaCirilloAC/paletter) by Andrea Cirillo does exactly this. However, the approach presented here is different and arguably much faster (with some caveats). 

The approach works as follows:

- First, the RGB cube (made of 255^3 colors) is transformed to LAB and then DIN99 color spaces. These spaces are perceptually uniform, which means that the Euclidean distance between colors in these spaces reflects better how humans perceive different colors. 

-  The DIN99 space is then summarized using k-means clustering, setting k = 2000. This means that every color with an RGB specification can be assigned to a cluster, and therefore to a cluster centroid; rounding the coordinates of a cluster centroid gives a mean color for each cluster. This summarized space and its clustering results are kept constant as they only need to be created once, and are distributed in the package.

- Colors are read from a picture and assigned to a summary color, quantifying the number of pixels that have been assigned to each color. 

- Colors are sorted, keeping the top 500 most abundant ones

- Optionally, colors can be filtered on high/low luminance and/or chroma (as defined by the HCL color space). Luminance is normally distributed and is filtered using mean +/- standard deviation(s), whereas chroma is filtered using median +/- median absolute deviation(s). 


- Sorted colors are projected again to the DIN99 space and k-means clustering is run again, setting a user-defined number of centers k, which corresponds to the palette size

- The most abundant color in the painting from each cluster is retained, so that every part of the color space is well represented.

<img width="458" alt="Screenshot 2022-06-04 at 6 19 48 PM" src="https://user-images.githubusercontent.com/21171362/171995101-3f5a11d8-06fe-4671-a11e-a18049a8583e.png">



The swatch will then have **k * m** colors, where **k** is the number of *clusters* (most diverse colors) and **m** is the number of *most representative colors per cluster* (hues within the cluster). This means that looking for 20 colors with k = 20 and m = 1 will yield a different result than k = 10 and m = 2, although the final number of colors will be the same; k = 20 gives more separated colors (as it will effectively partition the DIN99 space in 20 separate clusters), whereas k = 10 and m = 2 will select color couples from each cluster. The difference may be subtle in many cases but it allows users to fine tune the process.


<img width="1115" alt="Screenshot 2022-06-04 at 4 59 58 PM" src="https://user-images.githubusercontent.com/21171362/171992426-28209a0b-319d-452e-89a7-d2ef9a34537d.png">
<p align = "center">
Left: k = 20, m = 1; right: k = 10, m = 2
</p>

This procedure is quite fast mainly because of two tricks:

- the use of a pre-clustered reduced color space greatly reduces the feature space and speeds up the lookup 
- there is no need to optimize the palette by randomization, as the projection in DIN99 space followed by clustering affords a better color separation

The other positive aspect is that if you are using `swatcher` to analyze paintings _en masse_, they will all be assigned to a shared feature space of 1000 colors, making other types of analyses - such as embedding in a reduced space - feasible. Identifying swatches using `swatcher` also amounts to some sort of "feature selection" for statistical analyses of color in paintings. 


# Usage

Using `swatcheR` is easy. The two main functions are `analyzePictureCol()` to get the summary color distribution from a picture, and `getPalette()` to generate a palette from the distribution. A third function, `plotWithPal()` is a convenience function that wraps the two main functions and shows the picture and the resulting palette side by side. 

```
library(swatcheR)

analysis = analyzePictureCol(link = "https://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg/2560px-Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg")

head(analysis)
#550704 #170B04 #761B04 #1B1C04 #0E2404 #5D5F04 
     17   19676      47    8096    2476     698 
    
pal = getPalette(analysis, n = 20, sub = 1)
pal
 [1] "#BCC37C" "#B7BB9F" "#BCA834" "#94A9A4" "#BD9D23" "#8799B5" "#8C896F" "#758293" "#5B7593" "#5F7678" "#4A6593" "#596B2E" "#365398" "#344876"
[15] "#4C4519" "#414137" "#2E413B" "#443527" "#2F3539" "#1C2A5C"

plotWithPal(link = "https://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg/2560px-Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg", bg = "white")
```

<img width="596" alt="Screenshot 2022-06-04 at 5 55 41 PM" src="https://user-images.githubusercontent.com/21171362/171994345-3cfcf103-d3cd-4ff4-bae4-7696bcb33a4b.png">

As you can see the palette changes slightly, because the k-means search has a stochastic component to it, although this is mitigated by using 100 restarts.

The user can optionally decide whether to apply filters to color brightness (luminance) and saturation (chroma), and how much these filters remove in proportion to the luminance/chroma value distribution. 

## Different reference spaces

The package comes with a reference space built on 2000 clusters in DIN99 space. However, if you want to use a different reference space, you can build your own using `makeReferenceSpace()`:

```
lab_1500 = makeReferenceSpace(k = 1500, space = "LAB")
```

The available options at the moment are CIELab ("LAB") and DIN99 ("DIN99"). 

# Benchmarking

Coming soon...

# Acknowledgements
- Aditi Rajagopal for the package sticker
- Andrea Cirillo for the `paletter` package
- Ross Ihaka, Paul Murrell, Kurt Hornik, Jason C. Fisher, Reto Stauffer, Claus O. Wilke, Claire D. McWhite, Achim Zeileis for the `colorspace` package
- Jose Gama for the `colorscience` package
- Hadley Wickham, Winston Chang, Lionel Henry, Thomas Lin Pedesen, Kohske Takahashi, Claus Wilke, Kara Woo, and Hiroaki Yutani for the `ggplot2` package
- Thomas Lin Pedersen, Berendea Nicolae, Romain François for the `farver` package
- Simon Urbanek for the `png` and `jpeg` packages
- Hadley Wickham for the `reshape2` package



