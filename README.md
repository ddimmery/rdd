rdd 0.54
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

    Version: 0.54
    * Removed LICENSE to comport with CRAN policy.
    * Renamed VERSION to NEWS.
    
    Version: 0.53
    * Fixed bug in the estimation of LATEs away from zero.
    
    Version: 0.52
    * Fixed bug with sharp RD estimation.
	* Fixed bug when cutpoint not equal to zero.

    Version: 0.51
    * Fixed a bug which resulted in crashes of the summary.RD and plot.RD functions. 
	* Big changes in the handling of bandwidths. 
	* Now accepts a vector of bandwidths as inputs. 
	* Some outputs have been modified to reflect this change.
	
    Version: 0.50 
    * Added better error messages. 
	* Fixed bug in plot.RDD. 
	* Other minor code changes. 
	
    Version: 0.42
    * Submission issue. 
    
    Version: 0.41
    * Fixed submission issue. 
	
    Version: 0.4
    * Made compatible with CRAN. 
	
    Version: 0.31
    * Added plot.RD documentation. 
    
    Version: 0.3
    * Improved plot method, fixed some various bugs.
    
    Version: 0.22 
    * Fixed handling of missing values. 
    
    Version: 0.21 
    * Fixed bug in RDestimate when frame=TRUE and no covariates. 
    
    Version: 0.2 
    * Added separate function for kernel weights, S3 methods for plotting and summarizing RD objects.
	
    Version: 0.1 
    * First version. IKbandwidth, DCdensity and RDestimate created and working. 


