## FAQ

### Are there other ways to use DINEOF than with the compiled Fortran code?

If you are not needing to exploit all features of the compiled Fortran code, some wrappers or implementations in popular analysis software tools can be exploited:

* There is an R code available but not controlled, validated or maintained by our group: http://menugget.blogspot.be/2012/10/dineof-data-interpolating-empirical.html#uds-search-results
* See also DINEOF implemented in EOFimputation of metR : https://rdrr.io/github/eliocamp/metR/man/ImputeEOF.html
* Then there is an online version for use with standard satellite data: http://www.dineof.net/
* A matlab wrapper can be found here https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/+dineof/
* There also seems to be an ArcGIS implementation: http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1897232
* Finally an IDL wrapper was probably prepared : https://groups.google.com/forum/#!topic/dineof/PbkQwu8t7Z4


### Does the reconstruction take into account the time proximity of the images, i.e., does yesterday's image influence more than an image 6 months ago?

All the images have the same weight in the reconstruction. But in the recent versions, a filter is applied on the time-covariance matrix. This filtering provides a better consistence between successive images. Also note that even if there is no explicit information on time closeness, if the process is more correlated to recent events, the EOF decomposition will detect it and provide naturally a reconstruction which uses recent information rather than older ones.

### How many days/months/years of data are necessary for a reconstruction?

As a general rule, the longest time series of images will provide the best results, as the principal modes of variability will be well reproduced. However, if the memory usage is limited, then the time series can be shortened.

### Is it possible to reconstruct images with almost 100% of missing data?

Yes, it is possible, but the result may not be reliable. In general, images containing less than 5% of data do not provide useful information and might affect the final result.

### How to deal with clouds that are not well detected by the algorithms?

There is a _quality control_ implemented in `DINEOF`: it provides a list of suspect pixel values according to:
- the proximity with detected clouds
- the distance to the local median value
- the spatial coherence with the reconstruction

### Why the reconstruction does partially remove gradients or variability?

Because the reconstruction with `DINEOF` is obtained as a truncated EOF decomposition. There are higher-order modes that are not considered for the final reconstruction, as the reconstruction with this number of modes does not provide the minimal cross-validation error.

It is also possible to keep the original pixels and only fill in the clouds. Doing so, the gradients and the variability is conserved but there might be artifacts at the edges of the filled regions.

### Why isn't there any graphical user interface yet ?

`DINEOF` is intended for users who need to repeat calculation a number of times (be it because they are scientists who will need to test different parameter values on reconstructions, or operational users aiming a recalculating regularly new analysis). In this case, rather than clicking every time across a series of menus, a command line approach with batch files is much more efficient.

### Can I use DINEOF with 1D or 3D fields instead of 2D images?

Yes, as `DINEOF` does not know anything about spatial position and even in the 2D case simply stacks all values of a given moment into a column vector. The same can be done for 1D and 3D very easily.
