

This code is for multivariate time-series segmentation 
based on the mixture of Probabilistic Principal Component Analysis.

The main program is 'simpcamerge'. You can use our synthetic dataset
or your own datasets. You must set the number of principal components
based on the screeplot.

The code of the Bottom-Up segmentation algorithm can be found in 
'pcaseg' and 'pcaresid'.

The 'ppcamod' is the code for clustering. The program called 'compat' is able
to compute the compatibility matrix based on its results. The merging of 
compatible clusters is evaluated by 'mergeclust'.
