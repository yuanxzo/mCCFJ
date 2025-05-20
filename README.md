# mCCFJ
A MATLAB Package for calculating seismic ambient noise cross-correlation and frequency-Bessel transformation.

> In this package, GPU acceleration is available.


# Usage

Quickly understand its usage by typing the following command in Matlab's command line window

```
>> mCCFJ.Help("mCCFJ")  % English
Help Document - mCCFJ
 
Program Version: 1.2.0
Configuration Requirements:
  1. Matlab R2022b or later
  2. Python with obspy
  3. Matlab R2022b or later can successfully call Python
Document Usage:
     methods(mCCFJ);
     mCCFJ.Help('FunctionName', 'En');
     mCCFJ.Help('FunctionName', 'Zh');
Release Date: 2025-05-20
        
All Rights Reserved. Help Document - End


>> mCCFJ.Help("mCCFJ","Zh")  % Chinese
帮助文档 - mCCFJ
 
程序版本：1.2.0
配置要求：
        1. Matlab R2022b或更新版本
        2. 配置有obspy的Python
        3. Matlab R2022b或更新版本能成功调用Python
文档用法:
          methods(mCCFJ);
          mCCFJ.Help('FunctionName', 'En');
          mCCFJ.Help('FunctionName', 'Zh');
发布日期：2025-05-20
        
All Rights Reserved. 帮助文档 - 结尾
```

This package contains the following functions, which are displayed with the command "methods(mCCFJ)"
```
>> methods(mCCFJ);
Static methods:

Help       correlate  distance   filtering  transform 
```
Among them, "Help" is the help document program, "correlate" is the cross-correlation calculation program, "distance" is the cross-correlation distance calculation program, "transform" is the dispersion analysis program, and "filtering" is the k-filter program for cross-correlation. Their detailed usage can be understood by typing "mCCFJ.Help('FunctionName');", such as "mCCFJ.Help('correlate');"
