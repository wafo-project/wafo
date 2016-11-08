function d = cc2lc(cc,def,plotflag,sa)
%CC2LC Calculates the number of upcrossings from a cycle count  
%
% CALL:  lc = cc2lc(cc,def,plotflag,sa);
%
%      lc = a two column matrix with levels and number of upcrossings.
%      cc = cycle count (possibly rainflow filtered).
%
%     def = 1, only upcrossings.
%           2, upcrossings and maxima (default).
%           3, upcrossings, minima, and maxima.
%           4, upcrossings and minima.
%
%plotflag = 0, no plotting.
%           1, plot the number of upcrossings overplotted
%              with Rice formula for the crossing intensity
%              for a Gaussian process (default).
%           
%
%     sa  = standard deviation of the process.
%           (Default estimates it from the number of upcrossings)
%
% Calculates the number of upcrossings from a cycle count, e.g.
% min2Max cycles or rainflow cycles.
%
% Example:
%   tp = dat2tp(load('sea.dat'));
%   mM = tp2mm(tp);
%   lc = cc2lc(mM);
%
%   close all;
%
% See also  plotlc, tp2lc

% NB! needs normpdf to be able to overplot Rice formula

% Tested on: matlab 5.3
% History:
% revised by PJ 09-Jan-2000
%   copy of mm2lc
% revised by pab 11.08.99
%   changed name from mm2cross to mm2lc
% revised by Per A. Brodtkorb 01.10.98
%   added: overplot the crossingspectrum with Rice formula for crossing  
%   intensity for a Gaussian process



if nargin<4
  sa=[]; % unknown stdev is default
end

if nargin<3||isempty(plotflag)
  plotflag=1; %default
end
if nargin<2||isempty(def)
 def=2; % default
end

if ((def<1) || (def>4))
  error('def must be one of (1,2,3,4).')
end

index=find(cc(:,1) <= cc(:,2));

if isempty(index)
  error('Error in input cc.')
end

cc=cc(index,:);
ncc=length(cc);

minima=[cc(:,1)  ones(ncc,1) zeros(ncc,1) ones(ncc,1)];
maxima=[cc(:,2) -ones(ncc,1) ones(ncc,1) zeros(ncc,1)];

extremes=[maxima; minima];
[temp index]=sort(extremes(:,1));
extremes=extremes(index,:);

ii=1;
n=length(extremes);
extr=zeros(n,4);
extr(1,:)=extremes(1,:);
for i=2:n
  if extremes(i,1)==extr(ii,1);
    extr(ii,2:4)=extr(ii,2:4)+extremes(i,2:4);
  else
    ii=ii+1;
    extr(ii,:)=extremes(i,:);
  end
end
[xx nx]=max(extr(:,1));

if def==4 % This are upcrossings + minima
  d=[extr(1:nx,1) cumsum(extr(1:nx,2))];
  d(nx,2)=d(nx-1,2);
end

if def==1 % This are only upcrossings
  d=[extr(1:nx,1) cumsum(extr(1:nx,2)) - extr(1:nx,4)];
end

if def==3 % This are upcrossings + minima + maxima
  d=[extr(1:nx,1) cumsum(extr(1:nx,2)) + extr(1:nx,3)];
end

if def==2 % This are upcrossings + maxima
  d=[extr(1:nx,1) cumsum(extr(1:nx,2)) + extr(1:nx,3)-extr(1:nx,4)];
end

%% Plots are made by plotlc

if plotflag  
  plotlc(d,2,0,sa);
end

