%ITMKURS_LAB4 Script to computer exercises 4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of Measured Loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Estimation from Time Signal

load switchingload
plot(x0(:,1),x0(:,2))

n = 32; param = [-1 1.2 n]; % Define discretization
u = levels(param);     % Discrete levels
delta = u(2)-u(1)      % Discretization step
TP = dat2tp(x0,delta); % Get turning points and rainflow filter

TP0 = dat2tp(x0);
plot(x0(:,1),x0(:,2),TP0(:,1),TP0(:,2),TP(:,1),TP(:,2))
length(x0), length(TP0), length(TP)

RFC0 = tp2rfc(TP0); subplot(1,2,1), ccplot(RFC0)
RFC = tp2rfc(TP);   subplot(1,2,2), ccplot(RFC)
beta = 6;                % Damage exponent
Dam0 = cc2dam(RFC0,beta) % Damage 
Dam = cc2dam(RFC,beta)   % Damage, after rainflow filter 
Dam/Dam0

dtp = dat2dtp(param,TP,delta); % Discretized turning points 
tp = [dtp(:,1) u(dtp(:,2))'];  % Discretized turning points 
T = length(dtp);               % Number of turning points
clf, plot(x0(:,1),x0(:,2),TP(:,1),TP(:,2),tp(:,1),tp(:,2))
v=axis; hold on,
plot([v(1:2)],[u(2:end)'-delta/2 u(2:end)'-delta/2],'k:')
hold off, axis([v(1:2) param(1:2)])

rfc = tp2rfc(tp);        % Get rainflow cycles
dam = cc2dam(rfc,beta)   % Damage, after discretization & rainflow filter 
dam/Dam0

tz = [2 2; 46.40 1; 114.5 2; 161 1; 225.1 2; 270 1; 337.5 2; ...
384.8 1; 433.2 2; 600 1];
[xxd,xd,z] = splitload(dtp,tz);
plot(xxd{1}(:,1),xxd{1}(:,2))
plot(xxd{2}(:,1),xxd{2}(:,2))
hmmplot(xd(:,2),z,xd(:,1),[1 2],'','',1)

% Estimation of the Subloads

dtp1 = dat2tp(xxd{1});
[mM1,Mm1] = tp2mm(dtp1);
F1 = dcc2cmat(mM1,n);
dtp2 = dat2tp(xxd{2});
[mM2,Mm2] = tp2mm(dtp2);
F2 = dcc2cmat(mM2,n);
cmatplot(u,u,{F1 F2})

[G1s,h1] = smoothcmat(F1);
G1 = smoothcmat(F1,1,1.0,0);
[G2s,h2] = smoothcmat(F2);
G2 = smoothcmat(F2,1,0.8,0);
cmatplot(u,u,{G1s G2s; G1 G2})
cmatplot(u,u,{F1 F2; G1 G2})

%Estimation of the Regime Process}

N1 = length(dtp1), N2 = length(dtp2)
N12 = 4; N21 = 4;
p1=N12/N1; p2=N21/N2;
P = [1-p1 p1; p2 1-p2]  % P-matrix
statP = mc2stat(P)      % Stationary distribution

GG = {G1 []; G2 []};
[Grfc,mu_rfc] = smctp2rfm(P,GG);
cocc(param,RFC,Grfc)
Frfc = dtp2rfm(dtp(:,2),n);
cmatplot(u,u,{Frfc Grfc*T/2})

beta = 6;
Dam0, Dam, dam
damG = cmat2dam(param,Grfc,beta)*T/2
damG/Dam0
beta = 4;
Dmat = cmat2dmat(param,Frfc,beta);
DmatG = cmat2dmat(param,Grfc,beta)*T/2;
cmatplot(u,u,{Dmat,DmatG},3)

[xsim,zsim] = smctpsim(P,GG,T);
figure(1), hmmplot(u(xd(:,2))',z,1:length(xd),[1 2],'','',1)
figure(2), hmmplot(u(xsim)',zsim,1:T,[1 2],'','',1)

% Decomposition of a Mixed Rainflow Matrix

n = 16; param = [-1 1.2 n]; % Define discretization
u = levels(param);     % Discrete levels
delta = u(2)-u(1)      % Discretization step
TP = dat2tp(x0,delta); % Get turning points and rainflow filter

dtp = dat2dtp(param,TP,delta); % Discretized turning points 
tp = [dtp(:,1) u(dtp(:,2))'];  % Discretized turning points 
T = length(dtp);               % Number of turning points
rfc = tp2rfc(tp);        % Get rainflow cycles
beta = 6;
dam = cc2dam(rfc,beta)   % Damage, after discretization & rainflow filter 
dam/Dam0

Frfc = dtp2rfm(dtp(:,2),n);  % Observed rainflow matrix

tz = [2 2; 46.40 1; 114.5 2; 161 1; 225.1 2; 270 1; 337.5 2; ...
384.8 1; 433.2 2; 600 1];
[xxd,xd,z] = splitload(dtp,tz);
hmmplot(xd(:,2),z,xd(:,1),[1 2],'','',1)
dtp1 = dat2tp(xxd{1});
[mM1,Mm1] = tp2mm(dtp1);
F1 = dcc2cmat(mM1,n);
G1 = smoothcmat(F1,1,0.8,0);
dtp2 = dat2tp(xxd{2});
[mM2,Mm2] = tp2mm(dtp2);
F2 = dcc2cmat(mM2,n);
G2 = smoothcmat(F2,1,0.6,0);

known1.F = {G1 []; G2 []};   % known min-max and max-min matrices
init1.P = P;                 % initial guess of P-matrix
warning off                  % Don't display warnings
[Fest1,Est1] = estsmctp(Frfc,'P','ML',known1,[],init1);
Est1.P          % Estimated P-matrix
mc2stat(Est1.P) % Estimated stationary distribution

known3.Ffun = 'f_funm';      % Function for calculating a submodel
known3.trModel2X = 'tr_m2x'; % transform from Model to X-vector
known3.trX2Model = 'tr_x2m'; % transform from X-vector to model
known3.param = param;
init3.P = P;       % initial guess of P-matrix
M1.x0=[0.1 0.1]; M1.s=0.15; M1.lam=2; % submodel 1
M2.x0=[0.5 0.7]; M2.s=0.1; M2.lam=1;   % submodel 2
init3.M = {M1 M2}; % initial guess of Models for min-max matrices

OPTIONS(2)  = 1e-1;    % the termination tolerance for x;
OPTIONS(3)  = 1e-1;    % the termination tolerance for F(x);
[Fest3,Est3] = estsmctp(Frfc,'P,CalcF','ML',known3,[],init3);
Est3.P          % Estimated P-matrix
mc2stat(Est3.P) % Estimated stationary distribution
Est3.M{:}       % Estimated parameters in models

beta = 3:0.2:8;
Dam0 = cmat2dam(param,Frfc,beta)/(T/2); % Damage from load signal
FrfcEst1 = smctp2rfm(Fest1.P,Fest1.F);
Dam1 = cmat2dam(param,FrfcEst1,beta);   % Damage, scenario 1
FrfcEst3 = smctp2rfm(Fest3.P,Fest3.F);
Dam3 = cmat2dam(param,FrfcEst3,beta);   % Damage, scenario 3
plot(beta,Dam0,'b',beta,Dam1,'r',beta,Dam3,'g')
plot(beta,Dam1./Dam0,'r',beta,Dam3./Dam0,'g')
