function [Frm,paramM,paramR,paramA] = cmat2rmcmat(F,param)
%CMAT2RMCMAT  Converts a cycle matrix from Range-Mean to min-max format.
%
% CALL:  Frm = cmat2rmcmat(F);
%       [Frm,paramM,paramR,paramA] = cmat2rmcmat(F,param);
% 
% Input:
%   F      = Cycle matrix in min-max format.   [nxn]
%   param  = [a b n] ; defines discretization of levels.
% Output: 
%   Frm    = Cycle matrix in range-mean format.   [nxn]
%   paramM = [aM bM n] ; defines discretization of mean-values.
%   paramR = [aR bR n] ; defines discretization of ranges.
%   paramA = [aA bA n] ; defines discretization of amplitudes.
%
% NB! Due to different discretizations for Frm and F,
%     the mean values will not be identical.
%
% See also  rmcmat2cmat, cc2cmat

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 18-May-2000
%   Updated help text.
% Created by PJ (Pär Johannesson) 12-Apr-2000

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

% Initiate matrices
n = length(F);
Frm = zeros(n,n);

% Convert

for i = 1:n
  for j = i:n
    if F(i,j) ~= 0
      r = j-i+1;
      m = floor((i+j)/2);
      Frm(r,m) = F(i,j);
    end
  end
end

if ni>1
  paramM = param;  % Mean values
  paramR(3) = n;   % Ranges
  paramR(1) = (paramM(2)-paramM(1))/(n-1)/2;
  paramR(2) = paramR(1)+(paramM(2)-paramM(1));
  paramA = [paramR(1:2)/2 paramR(3)];  % Amplitudes
end
