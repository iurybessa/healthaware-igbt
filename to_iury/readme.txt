%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                               Description of files                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%F
FROM THE PAPER IT IS POSSIBLE TO COMBINE SOME FEATURES INTO ONE FEATURE USING REGRESSION TOOLS WITH THE REFERENCE BELOW, WHICH USED ANN 
TO IMPLEMENT MULTIPLE EXTRACTED FEATURES. THUS IT THINK IT CAN BE DONE, IT IS SOME SORT OF FEATURE FUSION LIKE IN PCA BUT THE FEATURE IS
THE INPUT AND THE FUSION DONE BY THE MODELING PROCESS IN SOME SENSE.

***Ali, Jaouher Ben et al. “Accurate bearing remaining useful life prediction based on
Weibull distribution and artificial neural network”. In: Mechanical Systems and
Signal Processing 56 (2015), pp. 150–172**


The folder contains 4 files ** I didnt do the acculmulation because the result is a convex fuction sort of and i guess it is not desirable**.

From the 4 files
*************************
features_Smoothen_features..... ===> contain a cell with the set of features, the first cell is the time.
from the second cell.
[RMS,Mean,Std,Skewness,kurtosis,...
   peak2peak,CrestFactor,ShapeFactor,ImpulseFactor,MarginFacto,Energy]
*******************
A monotocity , trendability and prognosability check was undertaken but for some strange reasons all the features have the same value  w.r.t to the prognosability and trendability.
 Thus you can select which features to combine from the monotonicity, enclosed is a figure (figure_show_monotonicity.jpg)in images showing the monotonicity
selcted features with high rates can be:::

[RMS,Mean,Std,MarginFactor]; .i.e (cell2,cell3,cell4,cell10)
*****************
Features_trig.mat ===>  contains a cell with all the features consdering the precursor signalk with smoothening and scaling done with a trignometric fn as done in the paper.
Basically i followed all the steps in teh paer except accummulate.

Apart from the first cell which is the time the rest as as follows


[RMS,Mean,Std,Skewness,kurtosis,...
   peak2peak,CrestFactor,ShapeFactor,ImpulseFactor,MarginFacto,Energy]

The monotonicity check was undertaken and the results can be seen from the figure (monotonicity_after_trig_smothening.jpg) in images.
selcted features with high rates can be:::

[[RMS,Mean,std,peak2peak...
  shapefactor...
   Mmean.Energy']]; .i.e (cell2,cell3,cell4,cell6,cell8,cell11)


**********************************
pca_for_normal==> Is the pca considering  smoothened selected features based on monotonicity check.

***********************************
 PCA_trig ==> Is the pca considering the trignometric fuction and scaling undertaken on the precursor signals.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I did not undertake the cummulation as there all require a bais to give a good feature.

