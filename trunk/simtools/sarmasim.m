function [x,z,e]=sarmasim(C,A,m,s2,P,T,Tinit,x0,z0)
%SARMASIM   Simulates a switching ARMA-process.
%   The regime process z has the state space {1,2,...,r}.
%   The process is governed by the system equation
%     A(q,z(t)) * x(t) = C(q,z(t)) * (m(z(t)) + s(z(t)) * e(t) )
%   with  s = sqrt(s2).
%
% [x,z,e] = sarmasim(C,A,m,s2,P,T,Tinit,x0,z0);
%
% x     = Simulated switching ARMA-process.
% z     = Simulated regime process.
% e     = Innovation process.
%
% C     = Coefficients in C-polynomials.        [rxnc+1]
% A     = Coefficients in A-polynomials.        [rxna+1]
% m     = Means of subprocesses.                [rx1]
% s2    = Innovation variances.                 [rx1]
% P     = Transition matrix for regime process. [rxr]
% T     = Length of simulation.
% Tinit = Length of simulation. (Optional, default: Tinit=10*na)
% x0    = Initial state of process x. If not given,
%         it will start from zeroes.            [1xna]
% z0    = Initial state of regime process. If not given, it will start 
%         from the stationary distribution of the Markov chain.
%
% Example: Switching ARMA(4,2)-process (Example 5 in thesis)
%   p1=0.005; p2=0.003; P = [1-p1 p1; p2 1-p2];
%   C = [1.00 1.63 0.65; 1.00 0.05 -0.88];
%   A = [1.00 -0.55 0.07 -0.26 -0.02; 1.00 -2.06 1.64 -0.98 0.41];
%   m = [46.6; 7.4]*1e-3;
%   s2 = [0.5; 2.2]*1e-3;
%   [x,z]=sarmasim(C,A,m,s2,P,2000);
%   plothmm(x,z)

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

if nargin < 6
  error('More input arguments needed. See help.')
end
if nargin<7, Tinit = []; end
if nargin<8, x0 = []; end
if nargin<9, z0 = []; end

r=size(P,1);
ma=size(A,2)-1;
mc=size(C,2)-1;
m_max=max([ma mc]);

if ma < 1
  A =[ones(r,1) zeros(r,1)];
  ma = 1;
end
if mc < 0
  C = ones(r,1);
  mc = 0;
end

if isempty(Tinit)
  Tinit=10*ma;
end
T1 = Tinit+T;
z=mcsim(P,T1,z0);

mx = (sum(C,2)./sum(A,2).*m)';
x=zeros(T1,1);
if isempty(x0)
  x(1:ma) = ones(ma,1).*mx(z(1:ma))';
else
  x(1:ma) = x0;
end

e=randn(T1,1);
s=sqrt(s2);

for t=m_max+1:T1
  x(t)=-A(z(t),2:ma+1)*x(t-1:-1:t-ma) + C(z(t),:)*(m(z(t))+s(z(t))*e(t:-1:t-mc));
end

x = x(Tinit+1:T1);
z = z(Tinit+1:T1);


