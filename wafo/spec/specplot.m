function specplot(S,varargin)
%SPECPLOT Plot a spectral density 
%
% CALL:  specplot(S,plotflag,linetype) 
%
%    S       = an array of spectral density structs:
%              1D  (see dat2spec)
%              2D  (see createspec)
%   1D:
%   plotflag = 1 plots the density, S, (default)
%              2 plot 10log10(S)
%	       3 plots both the above plots 
%   2D:
%   Directional spectra: S(w,theta), S(f,theta)             
%   plotflag = 1 polar plot S (default)
%              2 plots spectral density and the directional 
%                spreading, int S(w,theta) dw or int S(f,theta) df
%              3 plots spectral density and the directional 
%                spreading, int S(w,theta)/S(w) dw or int S(f,theta)/S(f) df
%              4 mesh of S
%              5 mesh of S in polar coordinates  
%              6 contour plot of S
%              7 filled contour plot of S
%   Wavenumber spectra: S(k1,k2)
%   plotflag = 1 contour plot of S (default)
%              2 filled contour plot of S
%   lintype : specify color and lintype, see PLOT for possibilities.
%
% NOTE: - lintype may be given anywhere after S.
%
% Examples
%  S = demospec('dir'); S2 = mkdspec(jonswap,spreading);
%  specplot(S,2), hold on
%  specplot(S,3,'g')  % Same as previous fig. due to frequency independent spreading
%  specplot(S2,2,'r') % Not the same as previous figs. due to frequency dependent spreading
%  specplot(S2,3,'m')
%  % transform from angular frequency and radians to frequency and degrees
%  Sf = ttspec(S,'f','d'); clf
%  specplot(Sf,2),
%
%  close all
%
% See also  dat2spec, createspec, simpson

plotspec(S,varargin{:})