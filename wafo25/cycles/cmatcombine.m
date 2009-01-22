function [F,Lim,FF1,FF2] = cmatcombine(F1,F2,in3)

%CMATCOMBINE Combines two cycle matrices.
%
% CALL: F = cmatcombine(F1,F2,rangeLim);
%       [F,Lim,FF1,FF2] = cmatcombine(F1,F2,Lim);
%
%   F1       = Cycle matrix 1                             [n,n]
%   F2       = Cycle matrix 2                             [n,n]
%   rangeLim = Use F1 for ranges >= rangeLim.             [1,1]
%   Lim      = Limitations on where to use F1.   [struct array]
%     .range = Use F1 for ranges >= Lim.range             [1,1]
%     .min   = Use F1 for min <= Lim.min                  [1,1]
%     .max   = Use F1 for max >= Lim.max                  [1,1]
%
%   F        = Cycle matrix, combination of F1 and F2     [n,n]
%   FF1      = Cycle matrix 1, used part                  [n,n]
%   FF2      = Cycle matrix 2, used part                  [n,n]
%
% Combine the two cycle matrices, F1 and F2, into one matrix, F,
% according to the conditions given by rangeLim (or Lim).
%
% Example:
%   F1 = triu(ones(8),0)
%   F2 = 2*F1
%   [F,Lim,FF1,FF2]=cmatcombine(F1,F2,2)
%   Lim=[]; Lim.range=2; Lim.min=4; Lim.max=4;
%   [F,Lim,FF1,FF2]=cmatcombine(F1,F2,Lim)
%
% See also  cmatplot, cc2cmat

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Pär Johannesson) 24-Jul-2000
% Updated by PJ 29-Aug-2000

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(3,3,ni));

% Third argument - Limitation
if isnumeric(in3)
  Lim.range = in3;
else
  Lim = in3;
end

n=length(F1); % Size of matrices

% Treat Lim
if ~isfield(Lim,'range'), Lim.range = []; end
if ~isfield(Lim,'min'),   Lim.min = []; end
if ~isfield(Lim,'max'),   Lim.max = []; end

% Default values
if isempty(Lim.range), Lim.range = 0; end
if isempty(Lim.min),   Lim.min = n; end
if isempty(Lim.max),   Lim.max = 1; end


% Combine the rainflow matrices

J=meshgrid(1:n);  % Columns = maximum
I=J';             % Rows    = minimum
K1 = (J-I>=Lim.range) &(I<=Lim.min) & (J>=Lim.max); % Where to use F1?
FF1 = F1; FF1(~K1)=0; % Set FF1 to zero where F1 shall NOT be used
FF2 = F2; FF2(K1) =0; % Set FF2 to zero where F1 shall be used
F = FF1+FF2;          % Combination

% Old code
% Combine the rainflow matrices

%F=F2;
%for i = 1:n
%  for j = i:n
%    if (j-i>=Lim.range) &(i<=Lim.min) & (j>=Lim.max)
%      F(i,j) = F1(i,j);
%    end
%  end
%end
