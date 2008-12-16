function lc=tp2lc(tp,def,plotflag,sa)
%TP2LC  Calculates the number of upcrossings from the turning points.
%
% CALL: lc = tp2lc(TP,def,plotflag,sa);
%
%      lc = a two column matrix with levels and number of upcrossings. [mx2]
%      TP = the turning points.                       [nx2]
%
%     def = 1, only upcrossings.
%           2, upcrossings and maxima (default).
%           3, upcrossings, minima, and maxima.
%           4, upcrossings and minima.
%
%plotflag = 0, no plotting
%           1, plot the number of upcrossings overplotted
%              with Rice formula for the crossing intensity
%              for a Gaussian process (default).
%           
%
%     sa  = standard deviation of the process
%           (Default estimates it from the number of upcrossings)
%
% See also  plotlc

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Pär Johannesson) 09-Jan-2000

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,4,ni));

if ni<2, def=[]; end
if ni<3, plotflag=[]; end
if ni<4, sa=[]; end

% Get min-max cycles
mM = tp2mm(tp); 
% Get level crossings
%lc = cc2lc(mM,def,plotflag,sa);
lc = cc2lc(mM,def,0,sa);
