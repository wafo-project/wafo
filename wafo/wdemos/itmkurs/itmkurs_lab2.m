%ITMKURS_LAB2 Script to computer exercises 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lab2.m 
% Switching Markov Loads and Rainflow Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% Section: Markov Chains of Turning Points
%

% Model Definition

n = 32; param = [-1 1 n];                % Define discretization
u = levels(param);                       % Discrete levels
G = mktestmat(param,[-0.2 0.2],0.15,1);  % min-max matrix

% Simulation

T = 5000;                % Length of simulation (number of TP)
xD = mctpsim({G []},T);  % Simulate
x = u(xD)';              % Change scale to levels -1,..,1

t=1:200; plot(t,x(t))  

% Calculating the Rainflow Matrix

Grfc=mctp2rfm({G,[]});

subplot(1,2,1),cmatplot(u,u,G),axis('square')
subplot(1,2,2),cmatplot(u,u,Grfc),axis('square')

cmatplot(u,u,{G Grfc},3)

FrfcObs = dtp2rfm(xD,n);
cmatplot(u,u,{FrfcObs Grfc*T/2})

% Switching Markov Chains of Turning Points

% Model Definition

n=32; param = [-1 1 n];                   % Define discretization
u=levels(param);                          % Discrete levels
G1 = mktestmat(param,[-0.4 -0.3],0.15,1); % regime 1
G2 = mktestmat(param,[0.3 0.4],0.15,1);   % regime 2

p1=0.10; p2=0.05;
P=[1-p1 p1; p2 1-p2]  % Transition matrix
statP=mc2stat(P)      % Stationary distribution

% Simulation

T=5000;   % Length of simulation (number of TP)
[xD,z] = smctpsim(P,{G1 []; G2 []},T); % Simulate
x=u(xD)'; % Change scale to levels -1,..,1

t=1:400;
hmmplot(x(t),z(t),t,[1 2])         % Same colour
hmmplot(x(t),z(t),t,[1 2],'','',1) % Different colours

% Calculating the Rainflow Matrix

Grfc=smctp2rfm(P,{G1,[];G2,[]});
Grfc1=mctp2rfm({G1,[]});
Grfc2=mctp2rfm({G2,[]});
GrfcSum=statP(1)*Grfc1+statP(2)*Grfc2;

cmatplot(u,u,{Grfc1 Grfc2; GrfcSum Grfc})

% Model B (Optional)

G1B = mktestmat(param,[-0.1 -0.1],0.28,0.5);  % regime 1
G2B = mktestmat(param,[0 0],0.12,2);          % regime 2

[xDB,zB] = smctpsim(P,{G1B []; G2B []},T); % Simulate
xB=u(xDB)'; % Change scale to levels -1,..,1
hmmplot(xB(t),zB(t),t,[1 2],'','',1) % Different colours

GrfcB=smctp2rfm(P,{G1B,[];G2B,[]});
Grfc1B=mctp2rfm({G1B,[]});
Grfc2B=mctp2rfm({G2B,[]});
GrfcSumB=statP(1)*Grfc1B+statP(2)*Grfc2B;

% Observed Rainflow Matrix and Smoothing

TP = dat2tp([(1:T)' xD]);      % Turning points
RFC = tp2rfc(TP);              % Rainflow cycles
paramD = [1 n n];
FrfcObs0 = cc2cmat(paramD,RFC); % Observed rainflow matrix

FrfcObs = dtp2rfm(xD,n); % Observed rainflow matrix

cmatplot(u,u,{FrfcObs/(T/2) Grfc},1)

h=0.8; FrfcSmooth = smoothcmat(FrfcObs,1,h);
cmatplot(u,u,{FrfcObs/(T/2) FrfcSmooth/(T/2) Grfc},1)

% Level Crossings

mu = cmat2lc(param,Grfc);
muSum = cmat2lc(param,GrfcSum);
muObs = cmat2lc(param,FrfcObs/(T/2));
subplot(2,1,1), plot(mu(:,1),mu(:,2),muSum(:,1),muSum(:,2),'--')
subplot(2,1,2), plot(mu(:,1),mu(:,2),muObs(:,1),muObs(:,2),'--')

% Damage

beta = 4;
Dam = cmat2dam(param,Grfc,beta)
DamSum = cmat2dam(param,GrfcSum,beta)
DamObs = cc2dam(u(RFC),beta)/(T/2)

Dmat = cmat2dmat(param,Grfc,beta);
DmatSum = cmat2dmat(param,GrfcSum,beta);
subplot(1,2,1), cmatplot(u,u,Dmat), axis('square'), v=axis;
subplot(1,2,2), cmatplot(u,u,DmatSum), axis('square'), axis(v)

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

% Methods for Estimation (Optional)

known.F = {G1 []; G2 []};   % known min-max and max-min matrices
init.P = P;                 % initial guess of P-matrix
[Fest,Est] = estsmctp(Fobs,'P','ML',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','chi2',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','HD',known,[],init); Est.P
[Fest,Est] = estsmctp(Fobs,'P','KL',known,[],init); Est.P
