function [RFM,RFM1,res] = dtp2arfm(x,in2,in3,res0)
%DTP2ARFM Calculates asymmetric RFM from discrete turning points.
%
% CALL:            RFM = dtp2arfm(x,n)
%       [RFM,RFM1,res] = dtp2arfm(x,n,def)
%       [RFM,RFM1,res] = dtp2arfm(x,def,RFM0,res0)
%
% Output:
%   RFM   = Asymmetric rainflow matrix (residual included).      [n,n]
%   RFM1  = Asymmetric rainflow matrix (without resudual).       [n,n]
%   res   = Residual.                                [nres,1]/[nres,2]
%
% Input:
%   x     = Turning points (taking values 1,...,n). [T,1]/[T,2]
%   n     = Number of levels.
%   def   = Treatment of residual.
%           'up':   Count min-to-Max cycles,    (default)
%                   gives correct number of upcrossings.
%           'down': Count Max-to-min cycles, 
%                   gives correct number of downcrossings.
%           'CS':   Cloormann/Seeger method, 
%                   gives all closed hysterisis loops.
%                   This method is identical to the French AFNOR recommendation,
%                   and the ASTM standard (variant called simplified version).
%   RFM0  = Asymmetric rainflow matrix (without resudual).     [n,n]
%   res0  = Residual (taking values 1,...,n).       [nres0,1]/[nres0,2]
%
% Calculates the Asymmetric RainFlow Matrix (ARFM) for the sequence of 
% discrete turning points,  by using the so-called 4-point algorithm.
%
% It is possible to calculate ARFM for 'dtp' and continue the counting 
% from previously counted 'ARFM0' with residual 'res0'
%   [ARFM,ARFM1,res] = dtp2arfm(dtp,def,ARFM0,res0)
%
% Example:
%   x = load('sea.dat');                   % Load data
%   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2
%   ARFM = dtp2arfm(dtp,32);            % Calculate asymmetric rainflow matrix
%   cmatplot(u,u,ARFM,3); colorbar;      % Plot rainflow matrix
%
%   close all;
% 
% See also  dtp2arfm4p, dtp2rfm, dcc2cmat, tp2rfc, dat2tp

% Tested on Matlab 5.3
%
% History:
% Revised by PJ 26-Jul-2000
%   New input 'def'. 
%   Now supports ANFOR and ASTM standards for rainflow counting.
% Revised by PJ (P�r Johannesson) 12-Jan-2000
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1999

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(2,4,ni));

if ni == 2
  n = in2; % Number of levels
  def = [];
  RFM = zeros(n);
  res0 = [];
elseif ni == 3
  n = in2; % Number of levels
  def = in3;
  RFM = zeros(n);
  res0 = [];
else
  def = in2;
  RFM0 = in3;
  n = length(RFM0);
  RFM = zeros(n);
end

% Default value
if isempty(def), def_res='up'; end

% Calculate RFM and res
[RFM,res] = dtp2arfm4p(x,n,res0);

% Add previously counted cycles (if any)
if ni == 4
  RFM = RFM+RFM0;
end

% Two or more outputs ?
if no >= 2
  RFM1 = RFM;
end

% Treat residual
% Calculate RFM = RFM0 + 'cycles in res'
ARFC_res = res2arfc(res(:,end),def);
for i=1:size(ARFC_res,1)
  RFM(ARFC_res(i,1),ARFC_res(i,2)) = RFM(ARFC_res(i,1),ARFC_res(i,2)) + 1;
end

