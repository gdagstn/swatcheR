<img src="https://user-images.githubusercontent.com/21171362/171850532-ded8fd04-6233-4f68-b278-2a23bf1ee495.png" align="right" alt="" width="150" />

# swatcheR
Fast automatic selection of color swatches from images

# Installation

You can install `swatcher` through `devtools`:

```
devtools::install_github("gdagstn/swatcheR")
```

# Motivation and methodology
This package lets you automatically and quickly select a color palette (or **swatch**) from a picture in .jpg or .png format. It is based on hierarchical or k-means clustering, which means that there will be some variability in the results unless you set a seed for the random number generator. 

The idea is not new: the `paletter` [package](https://github.com/AndreaCirilloAC/paletter) by Andrea Cirillo does exactly this. However, the approach presented here is different and arguably much faster (with some caveats). 

The approach works as follows:

- First, the RGB cube (made of 255^3 colors) is transformed to LAB and then DIN99 color spaces. These spaces are perceptually uniform, which means that the Euclidean distance between colors in these spaces reflects better how humans perceive different colors. 

-  The DIN99 space is then summarized using k-means clustering, setting k = 4096 (square root of 256^3). This means that every color with an RGB specification can be assigned to a cluster, and therefore to a cluster centroid; rounding the coordinates of a cluster centroid gives a mean color for each cluster. This summarized space and its clustering results are kept constant as they only need to be created once, and are distributed in the package.

- Colors are read from a picture and assigned to a summary color, quantifying the number of pixels that have been assigned to each color. "Extreme" values, such as total white/total black and other colors for which the R, G, B values are either 255 or 0, are kept separately to accommodate for very sharp colors (e.g. a total white background).

- Colors are sorted, keeping the top 2000 most abundant ones by default - this can be changed by the user.

- Optionally, colors can be filtered on high/low luminance and/or chroma (as defined by the HCL color space). Luminance is normally distributed and is filtered using mean +/- standard deviation(s), whereas chroma is filtered using median +/- median absolute deviation(s). 

- Sorted colors are clustered in one of two possible ways:
     - **Hierachical clustering** (the default): colors are projected to CIELab space and their $\Delta E$ 2000 distance is measured. The distance matrix is used for hierachical clustering, and the dendrogram is cut at a height resulting in a user-defined number of clusters, corresponding to the palette size. 
     - **K-means clustering**: colors are projected again to the DIN99 space and k-means clustering is run again, setting a user-defined number of centers k, which corresponds to the palette size
     - an additional method, **classic K-means** (for lack of a better name) applies k-means clustering directly on the DIN99 space, bypassing summarization in the reduced color space. It is a slower method that can still yield good results, but does not allow to compare palettes on a shared space.

- The most abundant color in the painting from each cluster is retained, so that every part of the color space is well represented.

<img width="880" alt="Screenshot 2022-06-05 at 2 10 51 PM" src="https://user-images.githubusercontent.com/21171362/172037835-bb66539c-a905-4152-8a7a-0a62972ffcfc.png">

_San Matteo e l'Angelo (Saint Matthew and the Angel)_, Caravaggio, 1602

The swatch will then have **k * m** colors, where **k** is the number of *clusters* (most diverse colors) and **m** is the number of *most representative colors per cluster* (hues within the cluster). This means that looking for 20 colors with k = 20 and m = 1 will yield a different result than k = 10 and m = 2, although the final number of colors will be the same; k = 20 gives more separated colors (as it will effectively partition the DIN99 space in 20 separate clusters), whereas k = 10 and m = 2 will select color couples from each cluster. The difference may be subtle in many cases but it allows users to fine tune the process.

<img width="1147" alt="Screenshot 2022-06-05 at 2 24 24 PM" src="https://user-images.githubusercontent.com/21171362/172038264-710610a1-f9e3-45cf-9cd2-5be370a46a58.png">

_Wheat Field with Cypresses_, Vincent Van Gogh, 1889

This procedure is quite fast mainly because of two tricks:

- the use of a pre-clustered reduced color space greatly reduces the feature space and speeds up the lookup 
- there is no need to optimize the palette by randomization, as the projection in perceptual spaces space followed by clustering affords a better color separation

The other positive aspect is that if you are using `swatcher` to analyze paintings _en masse_, they will all be assigned to a shared feature space of 4096 colors, making other types of analyses - such as embedding in a reduced space - feasible. Identifying swatches using `swatcher` also amounts to some sort of "feature selection" for statistical analyses of color in paintings. 

## What about continuous palettes?

It is notoriously hard to interpolate colors that have a high variability in hue. For a given categorical palette several different continuous palettes can be automatically generated. `swatcheR` tries to find good continuous palettes borrowing some tricks from computational geometry and graph theory:

- First, a large (n > 100) palette is used as input, and projected in DIN99 space. 
- The palette is clustered via k-means (k set to the square root of the number of colors in the palette)
- a k-nearest neighbor graph is built on these clusters, using k = 2 neighbors
- only the largest component of the graph is retained, and a minimum spanning tree is calculated for the subgraph
- the longest path on the MST of the subgraph is retained
- a principal curve is fit in DIN99 space using only colors belonging to the clusters included in the longest MST path
- principal curve coordinates are converted to colors.

The clustering step introduces a high level of randomness, which in this case is desirable, as many different continuous swatches can be built from a categorical one. 

This is the result of 12 random swatches from Caravaggio's _I bari_:

<img width="500" alt="Screenshot 2022-06-06 at 1 37 48 AM" src="https://user-images.githubusercontent.com/21171362/172063105-f524933b-17ef-4159-9de6-e9ee54a7d234.png">

_I bari (Cardsharps)_, Caravaggio, 1594

<img width="300" alt="Screenshot 2022-06-06 at 1 36 34 AM" src="https://user-images.githubusercontent.com/21171362/172063070-cd6dc50f-b00c-4fbd-b2ee-939109215ad8.png">

As you can see the yellow tones tend to appear frequently, but all swatches are different. 

# Usage

Using `swatcheR` is easy. The two main functions are `analyzePictureCol()` to get the summary color distribution from a picture, and `getPalette()` to generate a palette from the distribution. A third function, `plotWithPal()` is a convenience function that wraps the two main functions and shows the picture and the resulting palette side by side. 

```
library(swatcheR)

analysis = analyzePictureCol(link = "https://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg/2560px-Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg")

head(analysis)
#550704 #170B04 #761B04 #1B1C04 #0E2404 #5D5F04 
     17   19676      47    8096    2476     698 
    
pal = getPalette(analysis = analysis, n = 20, sub = 1)
pal
 [1] "#BCC37C" "#B7BB9F" "#BCA834" "#94A9A4" "#BD9D23" "#8799B5" "#8C896F" "#758293" "#5B7593" "#5F7678" "#4A6593" "#596B2E" "#365398" "#344876"
[15] "#4C4519" "#414137" "#2E413B" "#443527" "#2F3539" "#1C2A5C"

plotWithPal(link = "https://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg/2560px-Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg", bg = "white")
```

<img width="517" alt="Screenshot 2022-06-05 at 2 29 01 PM" src="https://user-images.githubusercontent.com/21171362/172038396-b13347d7-d1a5-4f4e-9efd-02540dfa6115.png">

_The Starry Night_, Vincent Van Gogh, 1889

The palette changes slightly, because the k-means search has a stochastic component to it, although this is mitigated by using 100 restarts.

The user can optionally decide whether to apply filters to color brightness (luminance) and saturation (chroma), and how much these filters remove in proportion to the luminance/chroma value distribution. 

For continuous palettes `getContinuousPalette()` can be used specifying the number of resulting palette (`tries`). It is a good idea to first generate a large palette:

```
pal = getPalette(analysis = analysis, n = 50, sub = 3)

cont_pal = getContinuousPalette(pal, tries = 5)

```
<img width="300" alt="Screenshot 2022-06-06 at 1 42 00 AM" src="https://user-images.githubusercontent.com/21171362/172063266-9ec9fa88-54ad-465d-bd7c-e6531cf94956.png">

## Hierachical clustering or k-means? 

There are some use cases in which HC outperforms k-means, and vice versa. For instance, if your picture has few colors, HC is the best choice:

<img width="764" alt="Screenshot 2022-06-05 at 2 21 02 PM" src="https://user-images.githubusercontent.com/21171362/172038112-beddeba8-662e-47ab-9d5e-a633a3fa3574.png">

If there are many different colors, k-means can have a higher range to capture nuances:

<img width="1077" alt="Screenshot 2022-06-05 at 8 00 25 PM" src="https://user-images.githubusercontent.com/21171362/172049432-e384c45f-446a-43cb-8b9d-64f514d7b524.png">

In the picture above, HC fails to capture the more saturated yellow component of the sunset, whereas k-means has a better coverage of the range. 

There are many parameteres that can be tuned to generate a satistfying palette.

## Different reference spaces

<img width="503" alt="Screenshot 2022-06-05 at 3 06 44 PM" src="https://user-images.githubusercontent.com/21171362/172039666-dfd2e000-237f-432c-ad7c-e35c4babf9a2.png">

A 3D representation of summary colors and their coordinates in the `din_4096` reference space, included in the package.

The package comes with a reference space built on 4096 clusters in DIN99 space. However, if you want to use a different reference space, you can build your own using `makeReferenceSpace()`:

```
lab_1500 = makeReferenceSpace(k = 1500, space = "LAB")
```

The available options at the moment are CIELab ("LAB") and DIN99 ("DIN99"). 


# Acknowledgements
- Aditi Rajagopal for the package sticker
- Andrea Cirillo for the `paletter` package
- Ross Ihaka, Paul Murrell, Kurt Hornik, Jason C. Fisher, Reto Stauffer, Claus O. Wilke, Claire D. McWhite, Achim Zeileis for the `colorspace` package
- Jose Gama for the `colorscience` package
- Hadley Wickham, Winston Chang, Lionel Henry, Thomas Lin Pedesen, Kohske Takahashi, Claus Wilke, Kara Woo, and Hiroaki Yutani for the `ggplot2` package
- Thomas Lin Pedersen, Berendea Nicolae, Romain François for the `farver` package
- Simon Urbanek for the `png` and `jpeg` packages
- Hadley Wickham for the `reshape2` package
- the `igraph` package [authors](https://cran.r-project.org/web/packages/igraph/AUTHORS)
- Trevor Hastie Andreas Weingessel, Kurt Hornik, Henrik Bengtsson, and Robrecht Cannoodt for the `princurve` package



