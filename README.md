rdd 0.50
========

This package provides the tools to undertake estimation in Regression 
 Discontinuity Designs. Both sharp and fuzzy designs are supported. Estimation is 
 accomplished using local linear regression. A provided function will utilize 
 Imbens-Kalyanaraman optimal bandwidth calculation. A function is also included to 
 test the assumption of no-sorting effects.
 
 
License
-------

[Apache License](http://www.apache.org/licenses/LICENSE-2.0.html)
 
Maintainer
----------

Drew Dimmery <[drewd@nyu.edu](mailto:drewd@nyu.edu)>
 
Changelog
---------

    --------------- 
    Version: 0.1 
    ----- 
    Changes: First version. IKbandwidth, DCdensity and RDestimate created and working. 
    --------------- 
    --------------- 
    Version: 0.2 
    ----- 
    Changes: Added separate function for kernel weights, S3 methods for plotting and summarizing RD objects. 
    --------------- 
    --------------- 
    Version: 0.21 
    ----- 
    Changes: Fixed bug in RDestimate when frame=TRUE and no covariates. 
    --------------- 
    --------------- 
    Version: 0.22 
    ----- 
    Changes: Fixed handling of missing values. 
    --------------- 
    --------------- 
    Version: 0.3 
    ----- 
    Changes: Improved plot method, fixed some various bugs.
    --------------- 
    --------------- 
    Version: 0.31 
    ----- 
    Changes: Added plot.RD documentation. 
    --------------- 
    --------------- 
    Version: 0.4 
    ----- 
    Changes: Made compatible with CRAN. 
    --------------- 
    --------------- 
    Version: 0.41 
    ----- 
    Changes: Fixed submission issue. 
    --------------- 
    --------------- 
    Version: 0.42 
    ----- 
    Changes: Submission issue. 
    --------------- 
    --------------- 
    Version: 0.50 
    ----- 
    Changes: Added better error messages. Fixed bug in plot.RDD. Other minor code changes. 
    --------------- 