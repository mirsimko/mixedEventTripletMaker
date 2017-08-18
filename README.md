Mixed event code for the Lambda_c analysis
==========================================

Based on the D0 mixed event code written by Michael Lomnitz

Maintainer:
  Miroslav Simko (msimko@bnl.gov)

How to install
--------------

```bash
git clone https://github.com/mirsimko/auau200GeVRun14.git
```
This version of auau200GeVRun14 is needed. All the code written on the original
auau200GeVRun14 by Jochen Thaeder, Mustafa Mustafa, and Xin Dong, should be 
compatible with this version as well.
```bash
git clone https://github.com/mirsimko/mixedEventTripletMaker.git
mkdir StRoot
cd StRoot
ln -s ../mixedEventTripletMaker/StRoot/StPicoMixedEventTripletMaker
ln -s ../auau200GeVRun14/StPicoHFMaker
ln -s ../auau200GeVRun14/StPicoCutsBase
ln -s ../auau200GeVRun14/StPicoPrescales
ln -s ../auau200GeVRun14/StPicoKFVertexFitter
mkdir macros
cd macros
ln -s ../mixedEventTripletMaker/StRoot/runPicoMixedEventTriplets.C
ln -s ../auau200GeVRun14/macros/loadSharedHFLibraries.C
```

Add mixed event maker classes to the `loadSharedHFLibraries.C` macro 
in `auau200GeVRun14/StRoot/macros` by adding a line
```C++
gSystem->Load("StPicoMixedEventTripletMaker");
```

Correct versions of picoDstMaker and StRefMultCorr are needed as well.

