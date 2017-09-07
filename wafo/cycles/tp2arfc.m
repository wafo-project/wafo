function [ARFC,ARFC1,res,def] = tp2arfc(x,def,ARFC0,res0)
%TP2ARFC Calculates asymmetric rainflow cycles from turning points.
%
% CALL:  [ARFC,ARFC1,res] = tp2arfc(tp,def,ARFC0,res0);
%                    ARFC = tp2arfc(tp);
%
% Output:
%   ARFC    = Asymetric RFC (residual included).      [N,2]/[N,4]
%   ARFC1   = Asymetric RFC (without resudual).       [N1,2]/[N1,4]
%   res     = Residual.                               [nres,1]/[nres,2]
%
% Input:
%   tp       = Turning points.                         [T,1]/[T,2]
%   def      = Choice of definition of rainflow cycles   [struct array]
%   def.res  = Treatment of residual.
%              'up':   Count min-to-Max cycles,    (default)
%                      gives correct number of upcrossings.
%              'down': Count Max-to-min cycles, 
%                      gives correct number of downcrossings.
%              'CS':   Cloormann/Seeger method, 
%                      gives all closed hysterisis loops.
%                      This method is identical to the French AFNOR recommendation, 
%                      and the ASTM standard (variant called simplified version).
%   def.time = 0: Don't store time of max and min. (default)
%              1: Store the time when the maxima and minima occured.
%   ARFC0    = Asymetric RFC (without resudual).       [N0,2]/[N0,4]
%   res0     = Residual.                               [nres0,1]/[nres0,2]
%
% Calculates the asymmetric rainflow cycles (ARFC) for the sequence of 
% turning points,  by using the so-called 4-point algorithm.
%
% It is possible to split the signal into smaller parts, and calculate 
% ARFC part by part. It can be especially useful for long signals.
% We count the first part and for the second part we continue counting 
% from previously counted 'ARFC0' with residual 'res0':
%   [ARFC1,ARFC0,res0] = tp2arfc(tp(1:1000,:));      % Firts 1000 points
%   [ARFC2] = tp2arfc(tp(1001:end,:),[],ARFC0,res0); % Point 1001 to end
% This shall give the same result as (i.e. ARFC=ARFC2)
%   [ARFC] = tp2arfc(tp);                            % Calculate all at once
%   sum(ARFC~=ARFC2)                                 % Shall return  [0 0]
%
% This routine doesn't use MEX, Fortran or C codes, only matlab code.
%
% Example:
%   x = load('sea.dat'); tp=dat2tp(x);
%   ARFC=tp2arfc(tp);      % Default (min-to-Max cycles in residual)
%   ccplot(ARFC);
%
%   close all;
%
% See also  tp2arfc4p, findrfc, dtp2arfm, tp2rfc, tp2mm, dat2tp, rfcfilt

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 06-Jul-2005
%   Fixed error with def & mod to avoid warning i R14SP2.
% Revised by PJ 26-Jul-2000
%   New format of def.
% Revised by PJ 12-Jul-2000
%   Now calls 'tp2arfc4p' to calculate ARFC0 and res.
%   Input 'def_res'.
%   Now supports AFNOR and ASTM standards for rainflow counting.
% Revised by PJ 18-May-2000
%   updated help text.
% Revised by PJ 09-Jan-2000
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1999

% Check input arguments
ni = nargin;
no = nargout;
%error(nargchk(1,4,ni));
narginchk(1,4) 
[T,nn] = size(x);

if ni<2, def = [];   end
if ni<3, ARFC0 = []; end
if ni<4, res0 = [];  end

def0=def;
if ~isempty(def)
  if isnumeric(def)
    def=[]; def.time = def0;
  elseif ischar(def)
    def=[]; def.res = def0;
  elseif ~isstruct(def)
    def=[];
  end
end

% Set default values
if ~isfield(def,'res')
  def.res = 'up';
end
if ~isfield(def,'time')
  def.time = 0;
end

% Calculate ARFC0 and res
if def.time == 0
  [ARFC,res] = tp2arfc4p(x(:,1:nn),res0,def.time);
else
  if nn==1
    [ARFC,res] = tp2arfc4p([(1:T)' x(:)],res0,def.time);
  else
    [ARFC,res] = tp2arfc4p(x,res0,def.time);
  end
end

% Add previously counted cycles (if any)
if ~isempty(ARFC0)
  ARFC = [ARFC0; ARFC];
end

% Rainflow cycles without residual
if no>=2, ARFC1=ARFC; end

% Rainflow cycles plus cycles in residual
% ARFC = ARFC + 'cycles in res'
ARFC_res = res2arfc(res,def.res,def.time); % Treat residual
ARFC = [ARFC; ARFC_res];


