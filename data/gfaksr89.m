%GFAKSR89  Reconstructed surface elevation measured at Gullfaks C 24.12.1989. 
%   
%  CALL:  xn = load('gfaksr89.dat')
%  
% Size             :    39000 X 2
% Sampling Rate    :    2.5 Hz
% Device           :    EMI laser
% Source           :    STATOIL 
% Format           :    ascii, c1: time c2: surface elevation
% Description      :
%   This is a reconstructed version of the data in the GFAKS89.DAT file.
%   The following calls were made to reconstruct the data:
%
%                inds = findoutliers(gfaks89,.02,2,1.23);
%            gfaksr89 = reconstruct(gfaks89,inds,6);
%
%   The wave data was measured 24th December 1989 at the Gullfaks C platform  
%   in the North Sea from 17.00 to 21.20. The period from 20.00 to 20.20
%   is missing in the original data.  The  water depth of 218 m is  
%   regarded as deep water for the most important wave components.
%   There are two EMI laser sensors named 219 and 220. This data set is 
%   obtained from sensor 219, which is located in the Northwest
%   corner approximately two platform leg diameters away from   
%   the closest leg.  
%   Thus the wave elevation is not expected to be significantly 
%   affected by diffraction effects for incoming waves in the western sector.   
%   The wind direction for this period is from the south.
%   Some difficulties in calibration of the instruments have been reported
%   resulting in several consecutive measured values being equal or almost equal 
%   in the observed data set. 
% 
%   Hm0 = 6.8m, Tm02 = 8s, Tp = 10.5
% 
% See also  gfaks89 

% Revised pab 04.10.2000
%  - fixed some typing errors.
