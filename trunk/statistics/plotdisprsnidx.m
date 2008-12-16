function di=plotdisprsnidx(data,varargin)
%PLOTDISPRSNIDX Plot Dispersion Index vs thresholds
%
% CALL  DI = plotdisprsnidx(data,options)
%
%  DI     = Dispersion Index as function of u.
%  data = twocolumn matrix with sampled times in the first column
%               and values the second columns
% options = options structure defining estimation of DI. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 3).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),100))
%        .Tb   : Block period (same unit as the sampled times)  (default 1)
%        .alpha: Confidence coefficient (default 0.05)
%
%   PLOTDISPRSNIDX plot the Dispersion Index (DI) as function of threshold. 
%   The purpose of DI is to determine the threshold where the number of
%   exceedances in a fixed period (Tb) is consistent with a Poisson process.
%   For a Poisson process the DI is one. Thus the threshold should be so
%   high that DI is not significantly different from 1.
%
%   The Poisson hypothesis is not rejected if the estimated DI is between:
%
%   chi2(alpha/2, M-1)/(M-1)< DI < chi^2(1 - alpha/2, M-1 }/(M - 1)
%
%  where M is the total number of fixed periods -generally
%   the total number of years in the sample.
%
% Example
%  xn = load('sea.dat');
%  Ie = findpot(xn,0,5);
%  di = plotdisprsnidx(xn(Ie,:),'Tb', 100); % a threshold around 1 seems appropriate.
%
% See also disprsnidx, reslife, fitgenparrange

di = disprsnidx(data,varargin{:});

if ~isstruct(di)
  plot(di,'.')
end
