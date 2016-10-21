function h_norm = smoothcmat_hnorm(F,NOsubzero)

%SMOOTHCMAT_HNORM  Bandwidth selection for kernel smoothing of a cycle matrix. 
%
% CALL: h_norm = smoothcmat_hnorm(F);
%       h_norm = smoothcmat_hnorm(F,NOsubzero);
%
% Input:
% F       = Cycle matrix.           [nxn]
% NOsubzero=Number of subdiagonals that are zero
%           (Optional, Default = 0, only the diagonal is zero)
%
% Output:
% h_norm    = Selected bandwidth.
%
% This choice is optimal if the sample is from a normal distribution
% The normal bandwidth usualy oversmooths, therefore one should choose 
% a slightly smaller bandwidth, e.g.  h=0.7*h_norm
%
% See also  cc2cmat, tp2rfc, tp2mm, dat2tp.

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Pär Johannesson) 18-Oct-2000
%   from  'smoothcmat'

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni<2, NOsubzero=[]; end

if isempty(NOsubzero), NOsubzero=0; end

n = length(F);    % Size of matrix
N = sum(sum(F));  % Total number of cycles

d = 2;   % 2-dim
[I,J] = meshgrid(1:n,1:n);

% Choosing bandwidth
% This choice is optimal if the sample is from a normal distr.
% The normal bandwidth usualy oversmooths,
% therefore we choose a slightly smaller bandwidth

h0 = N^(-1/(d+4));
FF = F+F';
mean_F = sum(sum(FF).*(1:n))/N/2;
s2 = sum(sum(FF).*((1:n)-mean_F).^2)/N/2;
s = sqrt(s2);       % Mean of std in each direction
h_norm = s*h0;      % Optimal for Normal distr.
h = h_norm;         % Test
