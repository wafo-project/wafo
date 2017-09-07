function [RFM,res] = dtp2arfm4p(x,n,res0)
%DTP2ARFM4P Calculates asymmetric RFM from discrete turning points (4-point).
%
% CALL:  [ARFM,res] = dtp2arfm4p(dtp,n);
%        [ARFM,res] = dtp2arfm4p(dtp,n,res0);
%
% Output:
%   ARFM  = Asymmetric RFM (without residual).       [n,n]
%   res   = Residual.                               [nres,1]/[nres,1]
%
% Input:
%   dtp   = Turning points (taking values 1,...,n). [T,1]/[T,2]
%   n     = Number of levels.
%   res0  = Residual (taking values 1,...,n).       [nres0,1]/[nres0,1]
%
% Example:
%   x = load('sea.dat');                   % Load data
%   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2
%   [ARFM,res] = dtp2arfm4p(dtp,32);      % Calculate asymmetric rainflow matrix
%   cmatplot(u,u,ARFM,3);                  % Plot rainflow matrix
%   colorbar;  
%   % res 
%
%   close all;
%
% See also  dtp2arfm, dtp2rfm, dcc2cmat, tp2rfc4p, dat2tp

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 26-Jul-2000
%   Created from 'dtp2arfm'

% Check input arguments
ni = nargin;
no = nargout;
%error(nargchk(2,3,ni));
narginchk(2,3)
if ni < 3
  res0 = [];
end

[T,nn] = size(x);
RFM = zeros(n);

nres = length(res0);
res = zeros(2*n+1,nn);
if nres>0
  res(1:nres,:) = res0;
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
      RFM(res(nres-2,nn),res(nres-1,nn)) = RFM(res(nres-2,nn),res(nres-1,nn)) + 1;
      res(nres-2,1:nn) = res(nres,1:nn);
      nres = nres-2;
    else
      cycleFound = 0;
    end
  end
end

% Residual
res = res(1:nres,:);


