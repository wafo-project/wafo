function [ro_min,ro_max,QQ]=mctp2stat(Q)
%MCTP2STAT  Calculates the stationary distribution for a MCTP.
%
% CALL: [ro_min,ro_max] = mctp2stat(F);
%
% ro_min = Stationary distribution of minima.         [1xn]
% ro_max = Stationary distribution of maxima.         [1xn]
%
% F      = Cell array of min-max and max-min  
%          matrices matrices for MCTP.                {1x2}
% F{1,1} = min-Max matrix                             [nxn]
% F{1,2} = Max-min matrix                             [nxn]
%
% Examples: 
%    [G, Gh] = mktestmat([-1 1 32],[-0.2 0.2],0.15,1);
%    [ro_min,ro_max] = mctp2stat({G Gh});
%
% See also  

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 18-May-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1999

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,1,ni));

if isempty(Q{1,2})
  Q{1,2} = Q{1,1}';
end

% Stationary distribution (=ro) of local minima with transition matrix
% Qt = Q*Qh = "Transition matrix for min-to-min"

Qt = Q{1,1}*Q{1,2};
ro_min = mc2stat(Qt(1:end-1,1:end-1));  % Stationary distr., row vector  
ro_min = [ro_min 0];  % Minimum can't reach the highest level

% Stationary distribution (=roh) of local maxima with transition matrix
% Qt = Qh*Q = "Transition matrix for max-to-max"

Qth = Q{1,2}*Q{1,1};
ro_max = mc2stat(Qth(2:end,2:end));  % Stationary distr., row vector  
ro_max = [0 ro_max];  % Maximum can't reach the highest level
