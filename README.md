# Point Cloud Clustering: Constant size, Compact and Non-Overlapping (C2NO)

This is a MATLAB function that receives a Point Cloud, and outputs
a set of clusters with the same number of points, specified by the user.
If the total number of points is not multiple of the target number of points per cluster,
the user can specify if only 1 cluster should be incomplete with
the remaining points, or if all points should be equally distributed.
These incomplete clusters can be padded by repeating points.
The resulting clusters are compact and do not overlap each other.
A detailed description of the algorithm will be available soon.

# Requirements

* **MATLAB:** Implemented and tested in **R2018b** on Linux (Ubuntu 18.04).
For other operating systems some changes might be required (e.g., on Windows, use backslash instead of forward slash for paths).
In previous versions of MATLAB, some functions might not be available (e.g., isfolder() was introduced in R2017b, previous versions use isdir() instead).
* **Computer Vision System Toolbox**

# Instructions for use:
### Usage
* c2no(infile, outpcname, outfolder, ppc)
* c2no(infile, outpcname, outfolder, ppc, cnn)
* c2no(infile, outpcname, outfolder, ppc, cnn, nb)
* c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic)
* c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic, padding)
* [progress, clust_disp, clust_bb, sum_dists, iter_adj_rest] = c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic, padding, debug)

### Input arguments:
* *infile* - Input filename of the point cloud to cluster (PLY or PCD)
* *outpcname* - Output clusters basename
* *outfolder* - Output clusters directory path
* *ppc* - Target number of points per cluster

### Optional input arguments:
* *cnn* - Number of nearest neighboring clusters eligible for transfers
     and refinement (Default: 10)
* *nb* - Number of block divisions on each direction for Cluster
     Centroid Initialization step, which divides the point cloud
     into nb*nb*nb blocks (Default: 8)
* *mic* - Selects the approach for dealing with target cluster sizes
     when the total number of points in the point cloud is not
     multiple of ppc; TRUE uses the Minimally Incomplete Clusters
     case, where points are equally distributed through all
     clusters; FALSE uses the Single Incomplete Cluster case, where
     all but one cluster will have ppc points, and one cluster will
     have the remaining points (Default: TRUE)
* *padding* - Selects if padding is applied to incomplete clusters,
     when the total number of points in the point cloud is not
     multiple of ppc (Default: TRUE)
* *debug* - This option allows tracking the progress of the algorithm
     in the different stages, computing cluster dispersion metrics
     at each iteration, and saving the clusters at the end of each
     phase (Default: FALSE)

### Output arguments:
* *progress* - 2D Matrix with the number of points in each cluster for
     each iteration of the Iterative Clustering Phase
* *clust_disp* - 3D Matrix with the cluster dispersion (measured as the
     variance of cluster coordinates) for each of the 3 dimensions,
     for each cluster, at each iteration of the Iterative Clustering
     Phase (only if debug=TRUE)
* *clust_bb* - 3D Matrix with the cluster dispersion (measured as the
     range of cluster coordinates) for each of the 3 dimensions,
     for each cluster, at each iteration of the Iterative Clustering
     Phase (only if debug=TRUE)
* *sum_dists* - 1D Matrix with the sum of distances between points and
     respective centroids, averaged for all clusters, taken at each
     iteration of the Clustering Refinement Phase (only if
     debug=TRUE)
     iter_adj_rest - Indicates at which iteration of the Iterative
         Clustering Phase the adjacency restriction was lifted, or the
         value 0 if it was not (only if debug=TRUE)
