function recintro
% RECINTRO  Info about RECDEMO: A statistical procedure for reconstruction of 1D signals. 
%
% The files in this directory can reproduce the figures in the 
% article of 
% 
% Brodtkorb, P.A., Myrhaug, D. and Rue, H (1999)
% 'Joint distribution of wave height and wave crest velocity
%  from reconstructed data.'
%  in Proceedings of 9th ISOPE Conference, Vol III, pp 66-73
% 
% 
% The paper shows a method for reconstruction of missing and or spurious 
% data points of a recorded time series by using a transformed Gaussian
% process. The application of the result is illustrated by fitting a
% joint distribution to the maximum of the crest front velocity and crest back velocity
% (V = max(Vcf,Vcb)) and the corresponding wave height (H = Hd or Hu)
% obtained from the reconstructed data set. The joint distribution seem to be physically 
% sound and much better described using the reconstructed than the 
% original time series. (Note that the extracted data are normalized with 
% the rms values, i.e., v=V/Vrms and h=H/Hrms)
% The joint probability distribution given as the product of the 
% Rayleigh distribution of wave heights and the Weibull distribution 
% for the conditional crest velocity (v) given wave height (h) seem to be
% a reasonable model.
%
%  DATA: GFAKS89.DAT 
% ~~~~~~~~~~~~~~~~~~~
% The wave data considered was measured 24th December 1989 at the Gullfaks C platform  
% in the North Sea from 17.00 to 21.20. The period from 20.00 to 20.20
% is missing and contains NaNs.  The  water depth of  218 m is  
% regarded as deep water for the most  important wave components.
% There are two EMI laser sensors named 219 and 220. This data set is 
% obtained from sensor 219, which is  located in the Northwest
% corner  approximately two platform leg diameters away from   
% the closest leg.  
%   Thus the wave elevation is not expected  to be significantly 
% affected by  diffraction effects for incoming waves in the western sector.   
% The wind direction for this period is from the south.
% Some  difficulties in calibration of the instruments have been reported
% resulting in several consecutive measured values being equal or nearly equal 
% in the observed data set.
%  
% Hm0 = 6.8m, Tm02 = 8s, Tp = 10.5 Tp2=19
%
% Wind data from Statfjord A, 1989 are also given to justify the
% assumption that the data may be regarded as stationary.




% By pab 28.01.2000

%more on
help recintro
%more off
