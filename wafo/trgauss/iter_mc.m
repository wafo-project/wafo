function [F, rfc]=iter_mc(f0,f,k,epsilon)
% ITER_MC  Calculates a kernel of a MC given a rainflow matrix
%
%    Solves  f_rfc = f_xy + F_mc(f_xy) for f_xy.
%
%  Call: [fmM_k frfc_k]=iter_mc(frfc,fmM_0,k,eps)
%
%   fmM_k  = the solution to the equation frfc = fmM + F(fmM),
%   frfc_k = the rainflow matrix; frfc_k = fmM_k + F(fmM_k).
%            
%
%   frfc   = the rainflow matrix to be inverted,
%   fmM_0  = the first approximation to the Markov matrix, if not
%            specified  fmM_0=frfc,
%   k      = number of iterations, if not specified, k=1.
%   eps    = a convergence treshold, default value; eps=0.00001
%
% See also  iter, spec2cmat, mctp2rfm, mc2rfm

% References:
% Rychlik, I. (1996)
% 'Simulation of load sequences from Rainflow matrices: Markov method'
% Int. J. Fatigue, Vol 18, pp 429-438
%

% tested on matlab 5.2
% History:
% by ir 1995


if nargin < 2
   f=f0;
end
if nargin <3
   k=1;
end
if nargin <4
   epsilon=0.00001;
end
check0=1;
f0=fliplr(f0);
f=fliplr(f);

for i=1:k
  if check0
   f1=f;
   rfc=mc2rfc(f);
   f=f1+(f0-rfc);
   f=max(0,f);
   check0=sum(sum(abs(f1-f)))>epsilon;
   %check=[k-i+1, sum(sum(abs(f1-f)))]
 end
end

F=fliplr(f);
rfc=fliplr(mc2rfc(f));

