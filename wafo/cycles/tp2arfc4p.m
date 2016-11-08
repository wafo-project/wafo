function [ARFC,res] = tp2arfc4p(x,res0,def_time)
%TP2ARFC4P Calculates asymmetric rainflow cycles from turning points (4-point).
%
% CALL:  [ARFC,res] = tp2arfc4p(tp)
%        [ARFC,res] = tp2arfc4p(tp,res0,def_time)
%
% Output:
% ARFC  = Asymmetric RFC (without residual).       [N,2]/[N,4]
% res   = Residual.                               [nres,1]/[nres,2]
%
% Input:
% tp       = Turning points.                         [T,1]/[T,2]
% res0     = Residual.                               [nres0,1]/[nres0,1]
% def_time = 0: Don't store time of max and min. (default)
%            1: Store the time when the maxima and minima occured.
%
% Calculates the rainflow cycles for the sequence of turning points, 
% by using the so-called 4-point algorithm.
%
% Calculate ARFC for turning points, starting from old residual 'res0'
%   [ARFC,res] = tp2arfc4p(tp,res0)
%
% This routine doesn't use MEX, Fortran or C codes, only matlab code.
%
% Example:
%   x = load('sea.dat'); 
%   tp = dat2tp(x);
%   [ARFC, res]=tp2arfc4p(tp);  % Default (min-to-Max cycles in residual)
%   ccplot(ARFC);
%   % res
%
%   close all;
%
% See also  tp2arfc, findrfc, dtp2arfm4p, tp2rfc, dat2tp, rfcfilt

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 2000-Jul-12
% Revised by PJ 04-Apr-2001
% - Added input def_time
% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(1,3,ni));
 
if ni<2, res0 = []; end
if ni<3, def_time = []; end

% Set default values
if isempty(def_time), def_time=0; end

[T,nn] = size(x);
ARFC = zeros(floor(T/2),2);
N = 0;

if def_time == 0, nn0=nn; else nn0=1; end

res = zeros(max([100 length(res0)]),nn);
if isempty(res0)
  nres = 0;
else
  nres = length(res0);
  res(1:nres,1:nn) = res0;
end

% Calculate ARFC and res
for i = 1:T
  nres = nres+1;
  res(nres,1:nn) = x(i,1:nn);
  cycleFound = 1;
  while cycleFound==1 && nres>=4
    if res(nres-1,nn) < res(nres-2,nn)
      A = [res(nres-1,nn) res(nres-2,nn)];
    else
      A = [res(nres-2,nn) res(nres-1,nn)];
    end
    if res(nres,nn) < res(nres-3,nn)
      B = [res(nres,nn) res(nres-3,nn)];
    else
      B = [res(nres-3,nn) res(nres,nn)];
    end
    %A = sort([res(nres-1) res(nres-2)]);
    %B = sort([res(nres) res(nres-3)]);
    if A(1) >= B(1) && A(2) <= B(2)
      N = N+1;
      arfc = [res(nres-2,nn:-1:nn0); res(nres-1,nn:-1:nn0)];      
      ARFC(N,1:2*(nn-nn0+1)) = arfc(:)';      
      res(nres-2,1:nn) = res(nres,1:nn);
      nres = nres-2;
    else
      cycleFound = 0;
    end
  end
end

% Counted rainflow cycles
ARFC = ARFC(1:N,:);

% Residual
res = res(1:nres,:);

