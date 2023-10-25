# color_repel
Repel visually similar colors away for colorblind users in various plots

For work (single cell RNA-seq) I make and look at countless scatterplots. Though most packages attempt to be colorblind-aware, I always find results uninterpretable when over a handful of colors are used. scatterHatch, adding different hatch patterns to clusters on top of colors, helps, but perhaps a simpler solution can be found -- avoid using visually similar colors next to each other (ie. on a UMAP, neighboring clusters should never be orange and yellow).

0. extract colors from plot object
1. generate distance matrix of categories (clusters on 2D plot, group on other types of plots)
2. generate distance matrix of colors, after conversion to CIELab space (and possibly various colorblindness conversion functions)
3. find optimal assignments of color to above groups/clusters
4. recolor
