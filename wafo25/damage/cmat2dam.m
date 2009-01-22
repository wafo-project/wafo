function D = cmat2dmat(param,F,beta,K)
%CMAT2DAM Calculates the total Palmgren-Miner damage of a cycle matrix.
% 
% CALL:  D = cmat2dam(param,F,beta,K);
%
% Output:
%   D     = Damage.                                            [1xm]
% Input: 
%   param = Parameter vector, [a b n], defines discretization.
%   F     = Cycle matrix.                                      [nxn]
%   beta  = Beta-values, material parameter                    [1xm]
%   K     = K-value, material parameter (Optional, Default: 1) [1x1]
% 
% The damage is calculated as 
%     D(i) = sum ( K * S^beta(i) ),  S = (max-min)/2
% 
% Example:
%   param = [-1 1 32]; F = mktestmat(param);
%   bv = 3:8; D = cmat2dam(param,F,bv); plot(bv,D,'x-')
%
% See also  cmat2dmat

% Tested on Matlab 6.0
%
% History:
% Revised by PJ 03-Nov-1999
% -  updated for WAFO
% Created by PJ (Pär Johannesson) 1997
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

% Calculate damage

n = length(F);
amp = cmat2amp(param,F);  % Histrogram of ranges

m=length(beta); D=zeros(1,m);

for i=1:m
  D(i) = K*sum((amp(:,1).^beta(i)).*amp(:,2));
end


