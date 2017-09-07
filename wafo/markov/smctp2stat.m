function [ro,ro_min,ro_max,Ro_min,Ro_max,QQ] = mctp2stat(P,F)
% SMCTP2STAT  Stationary distributions for a switching MCTP.
%
% CALL: [ro,ro_min,ro_max] = smctp2stat(P,F);
%       [ro,ro_min,ro_max,Ro_min,Ro_max] = smctp2stat(P,F);
%
% ro     = Stationary distribution of regime process.   [1xr]
% ro_min = Stationary distr. of minima for subloads. {r}[1xn]
% ro_max = Stationary distr. of maxima for subloads. {r}[1xn]
% Ro_min = Stationary distr. of minima for joint MCTP.  [1xnr]
% Ro_max = Stationary distr. of maxima for joint MCTP.  [1xnr]
%
% P      = Transition matrix for regime process.        [rxr]
% F      = Cell array of min-Max and Max-min matrices   {rx2}
% F{i,1} = min-Max matrix, process i                    [nxn]
% F{i,2} = Max-min matrix, process i                    [nxn]
%
% See also  mctp2stat, smctp2joint

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 18-May-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1999

% Check input arguments

ni = nargin;
no = nargout;
%error(nargchk(2,2,ni));
narginchk(2,2)
% Define 

r = length(P);   % Number of regime states
n = length(F{1,1});  % Number of levels

% Check that the rowsums of P are equal to 1

P = mat2tmat(P);

% Normalize the rowsums of F{1,1},...,F{r,1} to 1
%  ==>  QQ{1,1},...,QQ{r,1}

for i = 1:r
  QQ{i,1} = F{i,1};
  QQ{i,1} = mat2tmat(QQ{i,1},1);
end

% Normalize the rowsums of F{1,2},...,F{r,2} to 1
%  ==>  QQ{1,2},...,QQ{r,2}

for i = 1:r
  
  if isempty(F{i,2})        % Time-reversible
    QQ{i,2} = F{i,1}';
  else                   % F{i,2} is given
    QQ{i,2} = F{i,2}; 
  end
    
  QQ{i,2} = mat2tmat(QQ{i,2},-1);

end

% Stationary distribution (=ro) for regime process.

ro = mc2stat(P);

% Stationary distribution (=ro) of local minima with transition matrix
% Qt = Q*Qh = "Transition matrix for min-to-min"

for i = 1:r
  [ro_min{i},ro_max{i}] = mctp2stat(F(i,:));
end

[Q,QQ] = smctp2joint(P,F)

[Ro_min,Ro_max] = mctp2stat(Q);


