function [F,param] = rmcmat2cmat(Frm,paramM,paramR)
% RMCMAT2CMAT Converts cycle matrix from Range-Mean to min-max format.
%
% CALL:  F = rmcmat2cmat(Frm);
%       [F,param] = rmcmat2cmat(Frm,paramM,paramR);
% Output:
%   F      = Cycle matrix in min-max format.   [nxn]
%   param  = [a b n] ; defines discretization of levels.
% Input: 
%   Frm    = Cycle matrix in range-mean format.   [nxn]
%   paramM = [a b n] ; defines discretization of mean-values.
%   paramR = [a b n] ; defines discretization of ranges.
%
% NB! Due to different discretizations for Frm and F,
%     the mean values will not be identical.
%
% See also  cmat2rmcmat

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 18-May-2000
%   Updated help text.
% Created by PJ (Pär Johannesson) 12-Apr-2000

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,4,ni));

% Initiate matirces
n = length(Frm);
F = zeros(n,n);

% Convert

for r = 1:n
  for m = 1:n
    if Frm(r,m) ~= 0
      i = ceil(m-(r-1)/2);
      j = ceil(m+(r-1)/2);
      if i<1 | i>n, 
        warning(['minimum out of range (' num2str(i) '). Setting to boundary of [1,' num2str(n) '].']);
        i = min(max(i,1),n);
      end
      if j<1 | j>n, 
        warning(['Maximum out of range (' num2str(j) '). Setting to boundary [1,' num2str(n) '].']);
        j = min(max(j,1),n);
      end
      F(i,j) = Frm(r,m);
    end
  end
end

if ni>1
  param = paramM;
end
