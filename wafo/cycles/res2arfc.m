function [ARFC] = res2arfc(res,def,def_time)
%RES2ARFC Calculates asymmetric rainflow cycles for a residual.
%
% CALL:   ARFC_res = res2arfc(res);
%         ARFC_res = res2arfc(res,def,def_time);
%
% Output:
%   ARFC_res = Asymmetric RFC for residual.             [N,2]/[N,4]
%
% Input:
%   res      = Residual.                               [nres,1]/[nres,2]
%   def      = Treatment of residual.
%              'up':   Count min-to-Max cycles,    (default)
%                      gives correct number of upcrossings.
%              'down': Count Max-to-min cycles, 
%                      gives correct number of downcrossings.
%              'CS':   Cloormann/Seeger method, 
%                      gives all closed hysterisis loops.
%                      This method is identical to the French AFNOR recommendation, 
%                      and the ASTM standard (variant called simplified version).
% def_time   = 0: Don't store time of max and min. (default)
%              1: Store the time when the maxima and minima occured.
%
% Example:
%   x = load('sea.dat'); tp=dat2tp(x);
%   [ARFC,res]=tp2arfc4p(tp);      % Default (min-to-Max cycles in residual)
%   ARFC_res = res2arfc(res);      % Cycles in residual
%   plotcc(ARFC); hold on; plot(ARFC_res(:,1),ARFC_res(:,2),'r.'); hold off;
%
%   close all;
%
% See also  tp2arfc, tp2arfc4p, findrfc, dtp2arfm, tp2rfc, dat2tp, rfcfilt

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 25-Jul-2000
% Revised by PJ 09-Oct-2000
%   Some small corrections.
% Revised by PJ 04-Apr-2001
% - Added input def_time
% Correction by PJ 08-Nov-2001
%   Changed 
%     [ARFC,res_res] = tp2arfc4p(res2,def_time);
%   to
%     [ARFC,res_res] = tp2arfc4p(res2,[],def_time);
% Correction by PJ 13-Jun-2003
%   Exit funtion if less than two points in residual

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(1,3,ni));

if ni<2, def=[]; end
if ni<3, def_time = []; end

% Default value
if isempty(def), def='up'; end
if isempty(def_time), def_time=0; end

% Initiate
[nres,nn] = size(res);

% Exit funtion if less than two points in residual
if nres<2
    ARFC = [];
    return
end

% Count min-to-Max cycles, gives correct number of upcrossings
if strcmp(lower(def),'up')
  if res(2,nn)-res(1,nn)>0 % First point is a minimum?
    i_start=1;
  else
    i_start=2;
  end
    
  I = i_start:2:nres-1;
  if def_time == 0
    ARFC = [res(I,nn) res(I+1,nn)];
  else
    ARFC = [res(I,2:nn) res(I+1,2:nn) res(I,1) res(I+1,1)];
  end
  
% Count Max-to-min cycles, gives correct number of downcrossings
elseif strcmp(lower(def),'down')
  
  if res(2,nn)-res(1,nn)>0 % First point is a minimum?
    i_start=2;
  else
    i_start=1;
  end
  
  I = i_start:2:nres-1;
  if def_time == 0
    ARFC = [res(I,nn) res(I+1,nn)];
  else
    ARFC = [res(I,2:nn) res(I+1,2:nn) res(I,1) res(I+1,1)];
  end
  
% Cloormann/Seeger, gives all closed hysterisis loops  
elseif strcmp(lower(def),'cs')
  
  res2 = [res;res];           % Concatenate two residuals
  res2 = rfcfilter(res2,0,1); % Get turning points
  [ARFC,res_res] = tp2arfc4p(res2,[],def_time);
  
end

