function [lc,alpha] = mm2lc(mm,def,plotflag,sa)
%MM2LC Extracts level-crossing spectrum from min2Max cycles.   
%
%  CALL: [lc, alpha] = mm2lc(mM,def,plotflag,sa);
%
%      lc = two column matrix with levels and number of upcrossings,
%           i.e., level-crossing spectrum.
%   alpha = irregularity factor, approximately Tmaxima/Tz   
%      mM = min2Max cycles (possibly rainflow filtered).
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
%     sa  = standard deviation of the process
%           (Default estimates it from the number of upcrossings)
%
% Example: 
%   x = load('sea.dat'); tp = dat2tp(x); mM = tp2mm(tp);
%   lc = mm2lc(mM);
%
% See also  plotlc


% Tested on: matlab 5.3
% History:
% revised pab Feb2004  
%  - added alpha  
% revised pab 25.04.2001
% -speeded up the for loop even further.
% revised pab 19.04.2001
% - fixed a bug: a forgotten transpose.
% revised pab 30.12.2000
% - vectorized the for loop to speed up things
% revised by pab 11.08.99
% changed name from mm2cross to mm2lc
% revised by Per A. Brodtkorb 01.10.98
% added: overplot the crossingspectrum with Rice formula for crossing  
% intensity for a Gaussian process



%error(nargchk(1,4,nargin))
narginchk(1,4)
% Default values
if nargin<4,                    sa=[]; end      % unknown stdev is default
if nargin<3||isempty(plotflag),  plotflag=1; end % default plot final result
if nargin<2||isempty(def),       def=2; end      % default upcrossings && maxima 


if ((def<1) || (def>4))
  error('def must be one of (1,2,3,4).')
end

index = find(mm(:,1) <= mm(:,2));

if isempty(index)
  error('Error in input mM.')
end

cc     = mm(index,:); clear index
ncc    = length(cc);

minima = [cc(:,1)  ones(ncc,1) zeros(ncc,1) ones(ncc,1)];
maxima = [cc(:,2) -ones(ncc,1) ones(ncc,1) zeros(ncc,1)];

extremes   = [maxima; minima];
[tmp, ind] = sort(extremes(:,1));
extremes   = extremes(ind,:);

if 1,
  % pab 30.12.2000    
  % indices to matching entries.
  ind = (diff(tmp) == 0);
  % Create position mapping vectors
  tmp = [1;~ind];
  iy  = cumsum(tmp);
  ix  = [find(~ind);2*ncc];
  
  ind = find(ind).';   % added transpose (pab 19.04.2001)
  % Alternative call for finding ix:
  %ix = (1:2*ncc).';
  %ix(ind) = [];
else     
  %Alternatively, ix,iy and ind may be found by the following: (slower)
  [tmp, ix, iy]  = unique(extremes(:,1));  
  ind = find(diff(iy)==0)';
end

clear tmp
extr = extremes(ix,:);  % Keep only unique crossing levels
nx   = size(extr,1);
if any(ind)    
  if 1,
    % pab 25.04.2001 speeded up the for loop:
    l = diff([0; ix]);      % run lengths (ie, number of crossings)
    ind1 = find([1 diff(ind)>1]);
    for iz = ind1
      jy1 = ind(iz); 
      jx = iy(jy1);
      jy2 = jy1+l(jx)-2;
      extr(jx,2:4) = extr(jx,2:4) + sum(extremes(jy1:jy2,2:4),1);
    end
  else
    % Old call:  kept just in case (slow)
    for iz = ind
      extr(iy(iz),2:4) = extr(iy(iz),2:4) + extremes(iz,2:4);
    end
  end
end
clear extremes

switch def
  case 1,% Only upcrossings
    lc=[extr(1:nx,1) cumsum(extr(1:nx,2)) - extr(1:nx,4)];
  case 2,% Upcrossings + maxima
    lc=[extr(1:nx,1) cumsum(extr(1:nx,2)) + extr(1:nx,3)-extr(1:nx,4)];
  case 3,% Upcrossings + minima + maxima
    lc=[extr(1:nx,1) cumsum(extr(1:nx,2)) + extr(1:nx,3)];
  case 4,% Upcrossings + minima
    lc=[extr(1:nx,1) cumsum(extr(1:nx,2))];
    lc(nx,2)=lc(nx-1,2);
end

%% Plots are made by plotlc
if plotflag  
  plotlc(lc,2,0,sa);
end
if nargout>1
  cmax   = max(lc(:,2));
  alpha  = cmax/ncc;% approximately Tmaxima/Tz
end
return




 % Old call: slow  (kept just in case)
%   ii=1;
%   n=length(extremes);
%   extr=zeros(n,4);
%   extr(1,:)=extremes(1,:);
%   for i=2:n
%     if extremes(i,1)==extr(ii,1);
%       extr(ii,2:4) = extr(ii,2:4)+extremes(i,2:4);
%     else
%       ii=ii+1;
%       extr(ii,:) = extremes(i,:);
%     end
%   end
%   [xx nx]=max(extr(:,1));
