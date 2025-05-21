# mCCFJ
### A MATLAB Package for calculating seismic ambient noise cross-correlation and frequency-Bessel transformation.

In this package, these are two main functions, mCCFJ.correlate and mCCFJ.transform. 

*mCCFJ.correlate* is used to calculate the cross-correlation or cross-coherency function of seismic waveforms. To ensure the efficiency of calculation, the function is calculated in the frequency domain, and the cross-correlation between any two stations in the same time window is calculated in the way of matrix parallelism. If necessary, GPU acceleration can be used. 
      
*mCCFJ.transform* is used for dispersion analysis of cross-correlation function. We provide a variety of frequency-wavenumber domain transformation methods to deal with different data, which is up to you. Compared with previous methods, we provide an enhanced version of frequency-Bessel transform here, i.e. *spatial windowed frequency-Bessel transform*. This new method can make the energy of the dispersion spectrum more concentrated and reduce spatial artifacts, which is beneficial to the analysis of seismic wave phase velocity and attenuation. For details, refer to the following paper, and thank you for quoting it if the *mCCFJ* program brings convenience to your research.

      Yang, B., Meng, H., Yuan, S., & Chen, X. (2025). 
      Reliable Multimodal Attenuation Estimation of Surface Waves Using Diffuse Ambient Noise: Theory and Applications.
      Preprints. https://doi.org/10.22541/essoar.174120292.25931294/v1

This feature can be used by specifying ops.win='hamming_1', 'hamming_2' or 'hamming_half' in *mCCFJ.transform*. In order to ensure the efficiency of computing, GPU acceleration can be used when necessary.


# Usage

Quickly understand its usage by typing the following command in Matlab's command line window

```
>> mCCFJ.Help("mCCFJ")  % English
Help Document - mCCFJ
 
Program Version: 1.2.0
Configuration Requirements: Matlab R2022b or later versions are recommended
Document Usage:
     methods(mCCFJ);
     mCCFJ.Help('FunctionName', 'En');
     mCCFJ.Help('FunctionName', 'Zh');
        
Release Date: 2025-05-20
All Rights Reserved
        
Help Document - End


>> mCCFJ.Help("mCCFJ","Zh")  % Chinese
帮助文档 - mCCFJ
 
程序版本：1.2.0
配置要求：建议使用Matlab R2022b或更新版本
文档用法:
          methods(mCCFJ);
          mCCFJ.Help('FunctionName', 'En');
          mCCFJ.Help('FunctionName', 'Zh');
        
发布日期：2025-05-20
All Rights Reserved.
        
帮助文档 - 结尾
```

This package contains the following functions, which are displayed with the command "methods(mCCFJ)"
```
>> methods(mCCFJ);
Static methods:

Help  correlate  distances  filtering  transform 
```
Among them, "Help" is the help document program. You can know the details of other functions through "mCCFJ.Help('FunctionName')", such as "mCCFJ.Help('correlate')", so you can see the following
```
>> mCCFJ.Help("correlate")
Help Document - mCCFJ.correlate: Seismic waveform cross-correlation calculation program
 
Usage:
      CC=mCCFJ.correlate(wave, rloc, Fs);
      CC=mCCFJ.correlate(wave, rloc, Fs, options);
Required Input Parameters: 
      wave    At least two seismic waveform data, an array of size (npts,nsta), for example: wave=randn(npts,nsta), where npts is the number of waveform sampling points, nsta is the number of stations or waveform channels
      rloc    The location of the station, where rloc(1,1:nsta) is the longitude in latitude and longitude coordinates or the x-coordinate in rectangular coordinates,
                  rloc(2,1:nsta) is the latitude in latitude and longitude coordinates or the y-coordinate in rectangular coordinates
                  rloc(3,1:nsta) is the elevation or vertical coordinate in rectangular coordinates (can be ignored)
      Fs      Sampling rate of waveform data, unit Hz
Optional Input Parameters (options): 
      AX      The type of station location, 'latlon' means it is latitude and longitude coordinates, 'xyz' means it is rectangular coordinates, default 'latlon'
      CL      The length of the sliding window for correlation calculation, requires 0<CL<=npts, default CL=npts
      OL      The overlap length of the sliding window, 0<=OL<CL, default OL=0
      FM      The maximum frequency of interest, requires FM<=Fs/2, default FM<=Fs/2
      NT      Time domain normalization, 'No' means not executed, 'OneBit' means OneBit normalization, default 'No'
      NF      Frequency domain normalization, 'No' means not executed, 'PSD' means normalized by the average power spectral density of all stations in this time window, 'ABS', complete spectral whitening, means each station divided by its own modulus, default 'No'
      TP      The taper of Fourier transformation, ['No','Hann','tukeywin_5','tukeywin_10'], default 'tukeywin_5', i.e., 5% tukeywin window
      RR      The inter-station distance information calculated by the 'mCCFJ.distance' function, default [] means the program will automatically call the 'mCCFJ.distance' function to calculate
      DF      Downsample factor for the waveform spectrum, default 1 means no downsample. Downsampling can speed up the calculation, but it is necessary to pay attention to its impact on the time domain cross-correlation result
      GPU     Whether to use GPU to speed up the calculation during the calculation process, 'No' means not using, 'Yes' means using, default 'No'
Output Parameters (CC, is a structure): 
      CC.acf  Frequency domain autocorrelation function of each station data
      CC.ccf  Frequency domain cross-correlation function sorted by inter-station distance
      CC.cct  Time domain cross-correlation function sorted by inter-station distance
      CC.ccr  The first column is the inter-station distance, the second and third columns are the station retrieval corresponding to the inter-station distance, and the fourth column is the azimuth
      CC.freq Frequency sequence corresponding to the frequency domain cross-correlation function
      CC.time Time sequence corresponding to the time domain cross-correlation function
      CC.info Some additional information reserved during the calculation process, where CC.info.cstack records the effective stacking times of each column of CC.ccf
Help Document - End
```
