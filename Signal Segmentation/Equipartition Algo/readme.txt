 This code is a simple (not speed optimized) implementation of Simultaneous  Signal Segmentation 
and Modelling based on Equipartition Principle  [3] based on the papers [1-4]. 
You can use it only for non-commercial purposes. If you use it, please cite [3] or [1-4].
You can find more details in www.csd.uoc.gr/~cpanag and www.csd.uoc.gr/~tziritas 
Files: 
   runEPModelSelection.m: scipt that implements the method
   The runEPModelSelection.m calls the EP.exe.
   EP.exe can be created by compliling the .c files that are included in the same folder as .m files. 
   Otherwise, you can download EP.exe from https://dl.dropboxusercontent.com/u/11081708/ep.rar  

This is a general framework for simultaneous
segmentation and modelling of signals based on an
Equipartition Principle (EP). According to EP, the
signal is divided into segments with equal reconstruction
errors by selecting the most suitable model to
describe each segment. In addition, taking into account
change detection on signal model an efficient signal
reconstruction is also obtained. The model selection
concerns both the kind and the order of the model.
The proposed methodology is very flexible on different
error criteria and signal features.

[1] C. Panagiotakis, K. Athanassopoulos and G. Tziritas, The equipartition of curves,
 Computational Geometry: Theory and Applications, Vol. 42, No. 6-7, pp. 677-689, 2009.

[2] C. Panagiotakis and G. Tziritas, Signal Segmentation and Modelling based on Equipartition Principle, 
International Conference on Digital Signal Processing, 2009.

[3] C. Panagiotakis and G. Tziritas, Simultaneous Segmentation and Modelling of Signals based 
on an Equipartition Principle, 20th International Conference for Pattern Recognition  (ICPR), 2010.

[4] C. Panagiotakis, A. Doulamis and G. Tziritas, Equivalent Key Frames Selection Based on Iso-Content Principles, IEEE Trans. on Circuits and Systems for Video Technology, Vol. 19, No. 3, pp. 447 - 451, 2009. 