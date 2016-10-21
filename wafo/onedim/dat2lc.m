function lc = dat2lc(x,h,def)
%DAT2LC Extracts level-crossing spectrum from data, 
%       optionally rainflowfiltered. 
%
% CALL:  lc = dat2lc(x,h,def);
%
%  x  = two column data matrix with sampled times and values.
%  h  = a threshold; 
%       if  h<=0, then a sequence of turning points is used (default); 
%       if  h>0, then rainflow filtered cycles are used
% def = 1, only upcrossings.
%       2, upcrossings and maxima (default).
%       3, upcrossings, minima, and maxima.
%       4, upcrossings and minima.
%
%  lc = two column matrix with levels and number of upcrossings,
%        i.e., level-crossing spectrum.
%
% Example: 
%  x = load('sea.dat'); 
%  lc = dat2lc(x,0.2,1);
%  plotlc(lc)  
%  plot(lc(:,1),lc(:,2))
%
%  See also  dat2tp, mm2lc, dat2crossind

%
% Tested on: matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised jr 02.04.2001
%  - added example, updated info
% revised pab 30.12.2000
%  - added internal plotflag
% revised pab 24.11.2000
% by  Per A. Brodtkorb 11.08.99
%

error(nargchk(1,3,nargin))
plotflag=0;
if nargin<2||isempty(h),  h=0; end
if nargin<3||isempty(def),  def=2; end

tp=dat2tp(x,max(h,0));
[n m]= size(tp);
mM = [tp(1:2:n-1,m) tp(2:2:n,m)];

if plotflag
  lc = mm2lc(mM,def,1,std(x(:,2)));
else
  lc = mm2lc(mM,def,0);
end


