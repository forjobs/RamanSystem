# RamanSystem #

The **RamanSystem** is designed to process, analyze, and classify Raman spectrum. It contains numerous powerful algorithms that are executed using a simple graphical user interface.The source code is available so anyone can contribute to the development of the program.

![flow chart](https://github.com/forjobs/RamanSystem/blob/master/flow%20chart.jpg)

![graphical user interface](https://github.com/forjobs/RamanSystem/blob/master/graphical%20user%20interface%20of%20the%20Raman%20processing%20system.jpg)

## Requirements ##

The general system requirements are roughly the same as the MATLAB system requirements.


## Spectrum Files ##

The program uses simple text files to store Raman spectrum.

A spectrum file consists of a column of wavenumbers and a column of corresponding intensities separated by tab characters.

![selectdata interface](https://github.com/forjobs/RamanSystem/blob/master/selectdata%20interface.jpg)

Example data file is provided in the sample.txt.

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
  





