# RamanSystem #

The **RamanSystem** is designed to process, analyze, and classify Raman spectrum. It contains numerous powerful algorithms that are executed using a simple graphical user interface.The source code is available so anyone can contribute to the development of the program.

![flow chart](https://github.com/forjobs/RamanSystem/blob/master/flow%20chart.jpg)

[**The flow chart**](https://github.com/forjobs/RamanSystem/blob/master/flow%20chart.jpg "the flow of data through the system") displays the structure and data stream of system. The system is organized into four modules that execute specialized functions that operate upon data loaded in the system. The primary functions of the system include loading/initialization, preprocessing, recognition, and analysis/classification. The spectral data is imported from a special format file, and the corresponding spectral range and resolution of original spectrum signal are shown. After reading of data, some preprocessing methods are available, such as denosing, despiking, and baseline removal. If the information of sample class is unknown, analysis and classification is unique choice. The system provides three modeling methods, including PCA (principal Components Analysis), LDA (Linear Discriminant Analysis), and PLS (Partial Least Squares Regression).

![graphical user interface](https://github.com/forjobs/RamanSystem/blob/master/graphical%20user%20interface%20of%20the%20Raman%20processing%20system.jpg)

The system is developed in the form of [**graphical user interface **](https://github.com/forjobs/RamanSystem/blob/master/graphical%20user%20interface%20of%20the%20Raman%20processing%20system.jpg "graphical user interface")(GUI). The interface of the system is divided into several sections. In the top left area, the graph of the currently loaded original spectrum signal is shown. In the top right area, the preprocessed result of the spectrum signal is illustrated. In the bottom left corner, the plot of peak recognition is displayed. In the bottom right corner, there are four boxes in which the user is allowed to select and manipulate the spectral data to execute specialized functions as depicted in  [**the flow chart**](https://github.com/forjobs/RamanSystem/blob/master/flow%20chart.jpg "the flow of data through the system").

## Loading/Initialization module ##
![selectdata interface](https://github.com/forjobs/RamanSystem/blob/master/selectdata%20interface.jpg)

In the bottom right corner of the interface, there is a box in which the user is allowed to load and save the spectral data. In the “[selectdata](https://github.com/forjobs/RamanSystem/blob/master/selectdata%20interface.jpg)” window, the user can load the spectral data files, i.e. files with “.txt” extension. After a file has been loaded, the graph of currently selected original spectrum signal is illustrated in the top left area of the interface, and spectral range and resolution of original spectrum signal are shown in the box.

## Preprocessing module ##
The raw data is preprocessed by using Savitzky-Golay (SG) smoothing or wavelet denoising to remove spikes due to cosmic rays and reduce noise in each spectrum. It is also possible to employ SG smoothing followed by wavelet denoising or vice versa. The system employs three-point zero-order SG filter to reduce the noise and retain all important spectral bands. The SG filter is a moving window based on local polynomial fitting procedure. The size of moving window is three and polynomial order is zero in this method. Because the SG filter is a nonlinear weighted smoothing function, it is guaranteed that high frequency noise is well suppressed. Baseline correction is required to estimate and subtract the background fluorescence from each spectrum. A fully automated algorithm based on wavelet feature points and segment interpolation (AWFPSI) is applied to remove the baseline. At the end of the preprocessing stage, the results can be saved at the MATLAB workspace for use in the subsequent modules.
## Recognition module ##
A fully automated Raman peak recognition algorithm based on continuous wavelet transformation and local signal to noise ratio (CWTLSNR) is employed in this module. This algorithm extracts feature points through continuous wavelet transformation and recognizes peak through local signal to noise ratio. Through the wavelet transformation, the feature points in the original spectrum signal, such as peak, valley, starting point, ending point and inflection point, are classified. Combining the characteristics of wavelet feature points and the relationship between local signal to noise ratio (LSNR) and peak signal to noise ratio (PSNR), the position of the spectral peak can be recognized accurately. 
## Analysis/Classification module ##
The system provides three algorithms for analysis and classification as follows: PCA (principal Components Analysis), LDA (Linear Discriminant Analysis), and PLS (Partial Least Squares Regression). 
### Principle component analysis ###
To reduce redundant information and highlight important features, the system utilizes the method called PCA. PCA is a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components. Because only a few principal components are required to represent the majority of the total variance in the data set, PCA is an excellent way to reduce the dimensionality of the data.
### Linear Discriminant Analysis ###
Linear discriminant analysis (LDA) is a generalization of Fisher's linear discriminant, a method used in statistics, pattern recognition and machine learning to find a linear combination of features that characterizes or separates two or more classes of objects. The resulting combination may be used as a linear classifier or, more commonly, for dimensionality reduction before later classification. LDA attempts to find linear combinations (called discriminant functions) of the input variables that maximize between-group variation. Each function is given a discriminant score to determine how well it predicts group placement. Discriminant function scores are obtained by multiplying each discriminant function's coefficients with the corresponding points in the spectrum. For each spectrum, a classification score is computed for each group and the spectrum is assigned to the group with the highest score.
### Partial least squares analysis ###
Partial least squares (PLS) is a method for exploring patterns of covariation between two (and potentially more) blocks of variables. This method constructs new predictor variable (called components) as linear combinations of the original predictor variables. PLS is something of a cross between principal component analysis and multiple linear regression, therefore combining information about the variances of both the predictors and the responses, while also considering the correlations among them.