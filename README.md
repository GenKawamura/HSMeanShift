# MeanShift
R package to perform clustering of vector data and functional data using the mean shift algorithm.

### version 1.2.0
- implemented meanshift in Rcpp and improved its performance (approximately 30 times faster in msClustering method)

### version 1.1-1
- changes to the way "projectCurveWavelets" handles irregularly spaced data
- minor fixes to Vignette 2

### version 1.1-0
- added support for functional data via wavelet smoothing and thresholding via the "projectCurveWavelets" function
- included Vignette 1 - Clustering via the Mean Shift Algorithm
- included Vignette 2 - Clustering Functional Data via the Mean Shift Algorithm

### version 1.0-2
- fixed a bug that caused the functions "msClustering" and "bmsClustering" to produce an error if called with a value different than the default value for the argument "kernel"
- fixed a typo in the documentation of the function "msClustering": "options( mc.core=n.cores )" in the older documentation now correctly reads "options( mc.cores=n.cores )"
- other minor edits to the package documentation
