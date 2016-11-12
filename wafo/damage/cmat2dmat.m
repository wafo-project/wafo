function Dmat = cmat2dmat(param,F,beta,K)
%CMAT2DMAT Calculates (Palmgren-Miner) damage matrix from a cycle matrix.
%
% CALL:  Dmat = cmat2dmat(param,F,beta,K)
%
% Input: 
%   Dmat  = Damage matrix.                     [nxn]
% Output:
%   param = Parameter vector, [a b n], defines discretization.
%   F     = Cycle matrix.                      [nxn]
%   beta  = beta exponent.                     [1x1]
%   K     = K-value, material parameter (Optional, Default: 1) [1x1]
%
% Example:
%   param = [-1 1 32]; 
%   F = mktestmat(param);
%   Dmat = cmat2dmat(param,F,6);
%   cmatplot(Dmat);
%
%   close all;
%
% See also  cmat2dam, cmatplot, cc2cmat

% Tested on Matlab 6.0
%
% History:
% Revised by PJ  04-Jan-2000
% -  updated for WAFO
% Created by PJ (P�r Johannesson) 1997
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.0'


% Check input and otput

ni = nargin;
no = nargout;
error(nargchk(3,4,ni));

if ni < 4
  K=[];
end

% Set default values

if isempty(K)
  K = 1;
end

% Calculate damage matrix

n = length(F);
u=levels(param);
Dmat = zeros(n,n);

for i=1:n-1
  for j=i+1:n
    Dmat(i,j) = ((u(j)-u(i))/2)^beta*F(i,j);
  end
end


