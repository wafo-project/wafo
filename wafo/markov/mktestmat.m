function [F,Fh] = mktestmat(param,x0,s,lam,NOsubzero)
%MKTESTMAT   Makes test matrices for min-max (and max-min) matrices.
%
% CALL: [F,Fh] = mktestmat(param,x0,s,lam,NOsubzero)
%
% Input:
%   param  = Parameter vector, [a b n], defines discretization.
%   x0     = Center of ellipse. [min Max]             [1x2]
%   s      = Standard deviation. (0<s<infinity)       [1x1]
%   lam    = Orientation of ellipse. (0<lam<infinity) [1x1]
%            lam=1 gives circles.
%   NOsubzero = Number of subdiagonals that are set to zero
%               (-Inf: no subdiagonals that are set to zero)
%               (Optional, Default = 0, only the diagonal is zero)
%
% Output:
%   F      = min-max matrix.                          [nxn]
%   Fh     = max-min matrix.                          [nxn]
%
% Makes a Normal kernel (Iso-lines are ellipses).
% Each element in F =(F(i,j)) is
%   F(i,j) = exp(-1/2*(x-x0)*inv(S)*(x-x0)');
% where
%   S = 1/2*s^2*[lam^2+1 lam^2-1; lam^2-1 lam^2+1]
%
% The matrix Fh is obtained by assuming a time-reversible process.
% These matrices can be used for testing.
%
% Example:
%   [F,Fh] = mktestmat([-1 1 32],[-0.2 0.2], 0.25,1/2);
%   u = levels([-1 1 32]); 
%   cmatplot(u,u,F,3); axis('square');
%   [F,Fh] = mktestmat([-1 1 32],[-0.2 0.2], 0.25,1/2,-Inf);
%   cmatplot(u,u,F,3); axis('square');
%
%   close all;

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1997
%   Copyright (c) 1997 by P�r Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997


% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(0,5,ni));

if ni<1, param = []; end
if ni<2, x0 = []; end
if ni<3, s = []; end
if ni<4, lam = []; end
if ni<5, NOsubzero=[]; end

if isempty(param), param = [-1 1 32]; end
if isempty(x0), x0 = [1 1]*(param(2)+param(1))/2; end
if isempty(s), s = (param(2)-param(1))/4; end
if isempty(lam), lam = 1; end
if isempty(NOsubzero), NOsubzero=0;end

if isinf(NOsubzero), NOsubzero=-(param(3)+1); end
  
u = levels(param);
n = param(3);

% F - min-Max matrix

F=zeros(n,n);
S = 1/2*s^2*[lam^2+1 lam^2-1; lam^2-1 lam^2+1];

for i = 1:min(n-1-NOsubzero,n)
  for j=max(i+1+NOsubzero,1):n
    x = [u(i) u(j)];
    F(i,j) = exp(-1/2*(x-x0)*inv(S)*(x-x0)');
  end
end

% Fh - Max-min matrix
if no>1
  Fh = F';   % Time-reversible
end




