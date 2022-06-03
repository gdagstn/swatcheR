<img src="https://user-images.githubusercontent.com/21171362/171850532-ded8fd04-6233-4f68-b278-2a23bf1ee495.png" align="right" alt="" width="150" />

# swatcheR
Fast automatic selection of color palettes from images

# Motivation and methodology
This package lets you automatically and quickly select a color palette (or **swatch**) from a picture in .jpg or .png format. It is based on k-means clustering, which means that there will be some variability in the results unless you set a seed for the random number generator. 

The idea is not new: the `paletter`[package](https://github.com/AndreaCirilloAC/paletter) by Andrea Cirillo does exactly this. However, the approach presented here is different and arguably much faster (with some caveats). 

The approach works as follows:

- it reduces the RGB cube (made of 255^3 colors) to 1000 "summary colors" using k-means clustering, setting k = 1000. This means that every color with an RGB specification can be assigned to a cluster, and therefore to a cluster centroid; rounding the coordinates of a cluster centroid gives a mean color for each cluster. This reduced cube and its clustering results are kept constant as they only need to be created once, and are distributed in the package.

- then, it reads colors from a picture and assigns them to a summary color, quantifying the number of pixels that have been assigned to each color. 

- it sorts the colors keeping the top 500 most abundant ones

- optionally, it can filter on high/low luminance and/or chroma (as defined by the HCL color space). Luminance is normally distributed and is filtered using mean +/- standard deviation(s), whereas chroma is filtered using median +/- median absolute deviation(s). 


- projects these sorted colors to a perceptually uniform color space (DIN99), in which Euclidean distances correspond better to how humans perceive color differences

- runs k-means clustering again in DIN99 space, setting a user-defined number of centers k, which corresponds to the palette size

- it selects the most abundant color in the painting from each cluster, so that every part of the color space is well represented.

The swatch will then have **k * m** colors, where **k** is the number of *clusters* (most diverse colors) and **m** is the number of *most representative colors per cluster* (hues within the cluster). This means that looking for 20 colors with k = 20 and m = 1 will yield a different result than k = 10 and m = 2, although the final number of colors will be the same; k = 20 gives more separated colors (as it will effectively partition the DIN99 space in 20 separate clusters), whereas k = 10 and m = 2 will select color couples from each cluster. The difference may be subtle in most of the cases but it allows users to fine tune the process.


<p float="left">
  <img src="https://user-images.githubusercontent.com/21171362/171857049-cf10f09b-a1ee-4557-a0dd-4098671c1440.png" height="400" />
  <img src="https://user-images.githubusercontent.com/21171362/171857470-fa09a096-f1ea-4f55-bd66-5ddce60328c6.png" height="400" /> 
</p>


This procedure is quite fast mainly because of two tricks:

- the use of a pre-clustered reduced color space greatly reduces the feature space and speeds up the lookup 
- there is no need to optimize the palette by randomization, as the projection in DIN99 space followed by clustering affords a better color separation

The other positive aspect is that if you are using `swatcher` to analyze paintings _en masse_, they will all be assigned to a shared feature space of 1000 colors, making other types of analyses - such as embedding in a reduced space - feasible. Identifying swatches using `swatcher` also amounts to some sort of "feature selection" for statistical analyses of color in paintings. 

# Benchmarking





