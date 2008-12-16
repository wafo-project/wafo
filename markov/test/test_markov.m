%TEST_MARKOV Quick test of the routines in module 'markov'
%
% This script is used for a quick test of the routines in
% Modul 'markov':
%
% The following routines are tested:
%
% 
% MC - Markov chain
% - SIMMC - Simulates a Markov chain with state space {1,2,...,r}.
% - MC2RFM - Calculates the rainflow matrix/intensity for a Markov chain.
% - MC2STAT - Calculates the stationary distribution for a Markov chain.
%
% MCTP - Markov Chains of Turning Points
% - MCTPSIM - Simulates a Markov chain of turning points
% - MCTP2RFM - Calculates the rainflow matrix for a MCTP.
% - MCTP2STAT - Calculates the stationary distribution for a MCTP.
%
% SMC - Switching Markov Chains 
% - SMC2RFM - Calculates the rainflow matrix/intensity for a switching Markov chain.
%
% SMCTP - Switching Markov Chains of Turning Points
% - SMCTPSIM - Simulates a switching Markov chain of turning points,
% - PLOTHMM - plots a Hidden Markov Model.
% - SMCTP2RFM - Calculates the rainflow matrix for a SMCTP.
% - SMCTP2STAT - Stationary distributions for a switching MCTP.
%
% Decomposition of a Mixed rainflow matrix
% - ESTSMCTP - Estimate SMCTP model from an observed rainflow matrix.
%
% Simulation from Rainflow matrix
% - RFM2DTP  Reconstructs a sequence of turning points from a rainflow matrix.

% History:
% Created by PJ (Pär Johannesson) 18-May-2000
% Updated by PJ (Pär Johannesson) 07-Jul-2005

figure(1), clf

echo on

%MAT2TMAT - Converts a matrix to a transition matrix.

F = magic(5)
mat2tmat(F)
mat2tmat(F,1)
mat2tmat(F,-1)
mat2tmat(F,-1,-3)
mat2tmat(F,-1,0)

pause

%%%
% MC - Markov chain
%

% Model Definition - MC
n = 16; param = [-1 1 n];                % Define discretization
u = levels(param);                       % Discrete levels
[G,Gh] = mktestmat(param,[-0.2 0.2],0.2,1);  % min-max matrix
P = G+Gh;

%SIMMC - Simulates a Markov chain with state space {1,2,...,r}.
T = 5000;           % Length of simulation (number of TP)
xD = mcsim(P,T);  % Simulate
x = u(xD)';         % Change scale to levels -1,..,1

t=1:200; plot(t,x(t))  

pause

%MC2RFM - Calculates the rainflow matrix/intensity for a Markov chain.
Grfc=mc2rfm(P);
tp = rfcfilter(xD,0,1);
FrfcObs = dtp2rfm(tp,n);

cmatplot(u,u,{G Grfc FrfcObs/(T/2)},13)

pause

%MC2STAT - Calculates the stationary distribution for a Markov chain.
ro = mc2stat(P);

clf, plot(u,ro)

pause

%%%
% MCTP - Markov Chains of Turning Points
%

% Model Definition - MCTP
n = 16; param = [-1 1 n];                % Define discretization
u = levels(param);                       % Discrete levels
[G,Gh] = mktestmat(param,[-0.2 0.2],0.2,1);  % min-max matrix

%MCTPSIM - Simulates a Markov chain of turning points
T = 5000;                % Length of simulation (number of TP)
xD = mctpsim({G []},T);  % Simulate
x = u(xD)';              % Change scale to levels -1,..,1

t=1:200; plot(t,x(t))  

pause

%MCTP2RFM - Calculates the rainflow matrix for a MCTP.

Grfc=mctp2rfm({G,[]});
FrfcObs = dtp2rfm(xD,n);

cmatplot(u,u,{G Grfc FrfcObs/(T/2)},13)

pause

%MCTP2ARFM - Calculates the rainflow matrix for a MCTP.

Garfc=mctp2arfm({G,[]});
FarfcObs = dtp2arfm(xD,n);

cmatplot(u,u,{Garfc FarfcObs/(T/2)},13)

pause

%MCTP2STAT - Calculates the stationary distribution for a MCTP.
[ro_min,ro_max] = mctp2stat({G Gh});
[ro_min,ro_max] = mctp2stat({G []});

clf, plot(u,ro_min,u,ro_max)

pause

%%%%
% SMC - Switching Markov Chains 
% 

P=[0.9 0.1; 0.2 0.8]; 
param = [-1 1 32]; u = levels(param);
[F1,F2] = mktestmat(param,[-0.3 0.3],0.15,1,-Inf);

%SIMSMC - Simulates a Switching Markov chain with state space {1,2,...,r}.
[x,z] = smcsim(P,{F1 F2},400); 
plothmm(x,z,(1:length(z))',[1 2])         % Same colour

pause

%SMC2RFM - Calculates the rainflow matrix/intensity for a switching Markov chain.
Frfc = smc2rfm(P,{F1;F2});
cmatplot(u,u,Frfc)

pause

%%%%
% SMCTP - Switching Markov Chains of Turning Points
% 

% Model Definition

n=16; param = [-1 1 n];                   % Define discretization
u=levels(param);                          % Discrete levels
G1 = mktestmat(param,[-0.4 -0.3],0.15,1); % regime 1
G2 = mktestmat(param,[0.3 0.4],0.15,1);   % regime 2

p1=0.10; p2=0.05;
P=[1-p1 p1; p2 1-p2]  % Transition matrix
statP=mc2stat(P)      % Stationary distribution

% SMCTPSIM - Simulates a switching Markov chain of turning points,
T=5000;   % Length of simulation (number of TP)
[xD,z] = smctpsim(P,{G1 []; G2 []},T); % Simulate
x=u(xD)'; % Change scale to levels -1,..,1

%plothmm - plots a Hidden Markov Model.
t=1:400;
plothmm(x(t),z(t),t,[1 2])         % Same colour

pause

plothmm(x(t),z(t),t,[1 2],'','',1) % Different colours

pause

%SMCTP2RFM - Calculates the rainflow matrix for a SMCTP.
Grfc=smctp2rfm(P,{G1,[];G2,[]});
Grfc1=mctp2rfm({G1,[]});
Grfc2=mctp2rfm({G2,[]});
GrfcSum=statP(1)*Grfc1+statP(2)*Grfc2;

cmatplot(u,u,{Grfc1 Grfc2; GrfcSum Grfc})

pause

%SMCTP2STAT - Stationary distributions for a switching MCTP.
[ro,ro_min,ro_max] = smctp2stat(P,{G1,[];G2,[]});
[ro,ro_min,ro_max,Ro_min,Ro_max,QQ] = smctp2stat(P,{G1,[];G2,[]});

subplot(2,1,1), 
plot(u,ro_min{1},u,ro_max{1}), hold on
plot(u,ro_min{2},u,ro_max{2}), hold off
uu = u(floor(1:0.5:n+0.5));
subplot(2,1,2), plot(uu,Ro_min,uu,Ro_max)

pause

%%%
% Decomposition of a Mixed rainflow matrix
% Model and Estimation

n = 8; param = [-1 1 n];
M1.x0=[-0.4 -0.3]; M1.s=0.15; M1.lam=1; 
M2.x0=[0.3 0.4]; M2.s=0.15; M2.lam=1;
G1 = mktestmat(param,M1.x0,M1.s,M1.lam);
G2 = mktestmat(param,M2.x0,M2.s,M2.lam);
P=[1-0.1 0.1; 0.05 1-0.05]                % Transition matrix
[xD,z] = smctpsim(P,{G1 []; G2 []},5000); % Simulate
Fobs = dtp2rfm(xD,n);                     % observed mixed rainflow matrix

%ESTSMCTP - Estimate SMCTP model from an observed rainflow matrix.

known1.F = {G1 []; G2 []};   % known min-max and max-min matrices
init1.P = P;                 % initial guess of P-matrix
[Fest1,Est1] = estsmctp(Fobs,'P','ML',known1,[],init1);
Est1.P     % Estimated P-matrix

known3.Ffun = 'f_funm';      % Function for calculating a submodel
known3.trModel2X = 'tr_m2x'; % transform from Model to X-vector
known3.trX2Model = 'tr_x2m'; % transform from X-vector to model
known3.param =param;
init3.P = P;       % initial guess of P-matrix
init3.M = {M1 M2}; % initial guess of Models for min-max mat
[Fest3,Est3] = estsmctp(Fobs,'P,CalcF','ML',known3,[],init3);
Est3.P     % Estimated P-matrix
Est3.M{:}  % Estimated parameters in models

Pest_z = estmc(z,2)

mc2stat(Est1.P)
mc2stat(Est3.P)
mc2stat(Pest_z)
mc2stat(P)

pause

% Methods for Estimation (Optional)

known.F = {G1 []; G2 []};   % known min-max and max-min matrices
init.P = P;                 % initial guess of P-matrix
[Fest,Est] = estsmctp(Fobs,'P','ML',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','chi2',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','HD',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','KL',known,[],init); Est.P

%%%%%
% RFM2DTP  Reconstructs a sequence of turning points from a rainflow matrix.

x=load('sea.dat');
param = [-2 2 64]; n=param(3);
dtp0 = dat2dtp(param,x(:,2));
[RFM,RFM0,res0] = dtp2rfm(dtp0,n);
dtp = rfm2dtp(RFM0,res0);
plot(1:length(dtp0),dtp0,'b',1:length(dtp),dtp,'r')

echo off



