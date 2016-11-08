function [RFM,RFM1,res] = dtp2rfm(x,varargin)
%DTP2RFM Calculates rainflow matrix from discrete turning points.
%
% CALL:            RFM = dtp2rfm(dtp,n)
%       [RFM,RFM1,res] = dtp2rfm(dtp,n,def)
%       [RFM,RFM1,res] = dtp2rfm(dtp,def,RFM0,res0)
%
% Output:
%   RFM   = Rainflow matrix (residual included).    [n,n]
%   RFM1  = Rainflow matrix (without resudual).     [n,n]
%   res   = Residual.                               [nres,1]/[nres,2]
%
% Input:
%   dtp   = Turning points (taking values 1,...,n). [T,1]/[T,2]
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
%   RFM0  = Rainflow matrix (without resudual).     [n,n]
%   res0  = Residual (taking values 1,...,n).       [nres0,1]/[nres0,2]
%
% Calculates the rainflow matrix (RFM) for the sequence of discrete turning 
% points,  by using the so-called 4-point algorithm.
%
% It is possible to split the signal into smaller parts, and calculate 
% RFM part by part. It can be especially useful for long signals.
% We count the first part and for the second part we continue counting 
% from previously counted 'RFM0' with residual 'res0':
%   [RFM1,RFM0,res0] = dtp2rfm(dtp(1:1000,:),32);    % Firts 1000 points
%   [RFM2] = dtp2rfm(dtp(1001:end,:),[],RFM0,res0);  % Point 1001 to end
% This shall give the same result as (i.e. ARFM=ARFM2)
%   [RFM] = dtp2rfm(dtp,32);                         % Calculate all at once
%   sum(sum(RFM~=RFM2))                              % Shall return  0
%
% Example:
%   x = load('sea.dat');                   % Load data
%   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2
%   RFM = dtp2rfm(dtp,32);                 % Calculate rainflow matrix
%   cmatplot(u,u,RFM,3); colorbar;         % Plot rainflow matrix
%
%   close all;
% 
% See also  dtp2arfm, dcc2cmat, tp2rfc, dat2tp

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 26-Jul-2000
%   New input 'def'. 
%   Now supports AFNOR and ASTM standards for rainflow counting.
% Revised by PJ (P�r Johannesson) 12-Jan-2000
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1999

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(2,4,ni));

% Calculate asymetric RFM
if no < 2
  RFM = dtp2arfm(x,varargin{:});
else
  [RFM,RFM1,res] = dtp2arfm(x,varargin{:});
end

% Convert to symetric rainflow
RFM = RFM+RFM';
RFM = triu(RFM);
if no >= 2
  RFM1 = RFM1+RFM1';
  RFM1 = triu(RFM1);
end

