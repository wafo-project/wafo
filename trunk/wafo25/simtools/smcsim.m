function [x,z] = smcsim(P,Qc,T,init)
%SMCSIM  Simulates a Switching Markov chain with state space.
%
% CALL: [x,z] = smcsim(P,Q,T);
%       [x,z] = smcsim(P,Q,T,init);
%
% x    = Simulated switching Markov chain
% z    = Simulated Regime process
%
% P    = Transition matrix for regime process     [rxr]
% Q      = Cell array of transition matrices        {r,1}
% Q{i}   = Transition matrix for Markov chain i     [nxn]
% T    = Length of simulation.
% init.x0 = Initial state of process x. If not given, it will start from
%          the stationary distribution of minima given z(1).
% init.z0 = Initial state of regime process. If not given, it will start 
%          from the stationary distribution of the Markov chain.
%
% Simulates a Switching Markov chain with state space {1,2,...,n}. 
% The regime process has state space {1,2,...,r}.
%
% Example: Simulation of a switching Markov chain with two regime states.
%   P=[0.98 0.02; 0.05 0.95]; n=16; 
%   Q{1} = rand(n,n)*diag(exp(5*((n:-1:1)-1)/n)); 
%   Q{2} = rand(n,n)*diag(exp(5*((1:n)-1)/n)); % They will be normalized to row sum 1.
%   [x,z] = smcsim(P,Q,400); plothmm(x,z)
%   init.z0 = 2; init.x0 = [];
%   [x,z] = smcsim(P,Q,400,init); plothmm(x,z,[],[1 2],'','',1)
%   init.z0 = []; init.x0 = 6;
%   [x,z] = smcsim(P,Q,400,init); plothmm(x,z,[],[1 2],'','',1)

% %Example: Simulation of a Markov chain % Does not work
%   P=[0.9 0.1; 0.05 0.95];
%   x = smcsim(1,P,1000);
%   plot(x)

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 19-May-2000
%   updated for WAFO
%   Corrected method for simulating starting conditions.
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

% Check input arguments

ni = nargin;
%no = nargout;
error(nargchk(3,4,ni));

if ni < 4,  init = []; end

if isempty(init)
  init.x0 = [];
  init.z0 = [];
end

% Set constants
Zstr = '123456789';


r = length(P);     % Number of regime states
n = length(Qc{1}); % Number of levels

% Check that the rowsums of P are equal to 1

sumP = sum(P,2).';
if sum(sumP == 1) ~= length(P)
  warning('WAFO:SMCSIM','Rowsums of P not equal to 1. Renormalizing.');
  for i = 1:length(P)
    P(i,:) = P(i,:)/sumP(i);
  end
end

% Check that the rowsums of Qc{1},...,Qc{r} are equal to 1

for i = 1:r
  sumQi = sum(Qc{i},2).';
  if sum(sumQi == 1) ~= length(Qc{i})
    warning('WAFO:SMCSIM',['Rowsums of Q{' Zstr(i) '} not equal to 1. Renormalizing.']);
    for j = 1:length(Qc{i})
      Qc{i}(j,:) = Qc{i}(j,:)/sumQi(j);
    end
  end
end


% Make the transition matrix Q for the joint MC (X_k,Z_k)

Q = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  QQ = kron(Qc{z},P);
  Q(I+z,:) = QQ(I+z,:);
end


% Stationary distribution = ro of Q

ro = mc2stat(Q);

% Generate random numbers

e=rand(T,1);

% Start values

e0 = e(1);
if isempty(init.z0) && isempty(init.x0)
  x0z0 = min(find( e0<=cumsum(ro) ));
  x0 = floor((x0z0+1)/r);
  z0 = mod(x0z0-1,r)+1;
elseif isempty(init.x0)
  z0 = init.z0;
  rox0 = ro(z0:r:end); % Pick stat. distr. for regime z0
  rox0 = rox0/sum(rox0);
  x0 = min(find( e0<=cumsum(rox0) ));
elseif isempty(init.z0)
  x0 = init.x0;
  z0 = [];  % Start from stat. distr of P
else % Both z0 znd x0 are given
  x0 = init.x0;
  z0 = init.z0;
end

% Simulate Regime process

z = mcsim(P,T,z0);

% Simulate switching Markov chain

x=zeros(T,1);
x(1) = x0;   % First value

for k=2:T
  Pi = Qc{z(k)};
%  eval(['Pi = P' Zstr(z(k)) ';']);
  cumsumPi = cumsum(Pi,2);
  x(k) = min(find( e(k)<=cumsumPi(x(k-1),:) ));
end
