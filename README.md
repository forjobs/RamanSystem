# RamanSystem #

The **RamanSystem** is designed to process, analyze, and classify Raman spectrum. It contains numerous powerful algorithms that are executed using a simple graphical user interface.The source code is available so anyone can contribute to the development of the program.

![flow chart](https://github.com/forjobs/RamanSystem/blob/master/flow%20chart.jpg)

![graphical user interface](https://github.com/forjobs/RamanSystem/blob/master/graphical%20user%20interface%20.jpg)

## Requirements ##

The general system requirements are roughly the same as the MATLAB system requirements.


## Spectrum Files ##

The program uses simple text files to store Raman spectrum.

A spectrum file consists of a column of wavenumbers and a column of corresponding intensities separated by tab characters.

![selectdata interface](https://github.com/forjobs/RamanSystem/blob/master/selectdata%20interface.jpg)

Example data file is provided by the *.txt in the folder "data". The sample data for each function is provided in the corresponging subfolder, and should be added into the path when software runs.

## Features ##
Here are some of the features of the Raman Processing program:

### Preprocessing ###
- Savitzky-Golay (SG) smoothing or wavelet denoising
- automated algorithm based on wavelet feature points and segment interpolation (AWFPSI)

### Recognition ###
- automated Raman peak recognition algorithm based on continuous wavelet transformation and local signal to noise ratio(CWTLSNR)

### Analysis/Classification ###
- Principal component analysis
- Linear discriminant analysis
- Partial least squares analysis
  





