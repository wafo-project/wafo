function x=mcsim(P,T,x0)

%MCSIM   Simulates a Markov chain.
%
% CALL: x = mcsim(P,T)
%       x = mcsim(P,T,x0)
%
% x  = Simulated Markov chain.
%
% P  = Transition matrix.    [rxr]
% T  = Length of simulation.
% x0 = Initial state. (Default:  x0 will be from the stationary
%      distribution of the Markov chain.)
%
% Simulates a Markov chain with state space {1,2,...,r}
%
% Example:
%   P=[0.8 0.2; 0.1 0.9];
%   x=mcsim(P,100); plot(x)
%
% See also  mc2stat, smcsim

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 26-Jul-2000
%   updated help text.
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

if nargin < 3
  x0 = [];
end

% Check that the rowsums of P are equal to 1

sumP = sum(P,2).';
if sum(sumP == 1) ~= length(P)
  disp('Warning: Rowsums of P not equal to 1. Renormalizing.');
  for i = 1:length(P)
    P(i,:) = P(i,:)/sumP(i);
  end
end

x=ones(T,1);
e=rand(T,1);
cumsumP = cumsum(P,2);

% Initial state

if isempty(x0)    % From stationary distribution

  ro=mc2stat(P);

  x(1) = min(find( e(1)<=cumsum(ro) ));
%  x(1) = sum( cumsum(ro) <= e(1) ) + 1;

else

  x(1) = x0;   % Given

end

%---------------------

for k=2:T

  x(k) = min(find( e(k)<=cumsumP(x(k-1),:) ));

end
 
%  x(k) = sum( cumsumP(x(k-1),:) <= e(k) ) + 1;
