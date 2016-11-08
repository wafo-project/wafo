%TEST_CYCLES Quick test of the routines in module 'cycles'
%
%This script is used for a quick test of the routines in
% Modul 'cycles':
%
% The following routines are tested:
%
% Cycle counting
% - dat2tp
% - rfcfilter
% - tp2mm
% - tp2rfc
% - tp2arfc
% 
% Cycle matrix
% - cc2cmat
% - cmat2nt
% - nt2cmat
% - cmat2rmcmat
% - rmcmat2cmat
% - smoothcmat
% - tp2lc
% - cc2lc
% - cmat2lc
% - nt2lc
% - cc2amp
% - cmat2amp
%
% Discrete loads
% - dat2dtp
% - cc2dcc
% - dtp2rfm
% - dtp2arfm
% - dcc2cmat
%
% Extrapolation & Smoothing of RFM/CMAT/LC. 
% - rfmextrapolate
% - tpextrapolate
% - cmat2extralc  
% - fitgenpar_mld  
% - lc2rfmextreme
% - extralc     
% - smoothcmat   

% History:
% Created by PJ (P�r Johannesson) 18-May-2000
% Updated by PJ (P�r Johannesson) 07-Jul-2005

figure(1), clf

echo('on')

% Load test-signal: deep.dat
x = load([waforoot filesep 'wdemos' filesep 'itmkurs' filesep 'deep.dat']);

plot(x(:,1),x(:,2))

pause

%%%%
% Cycle counting
%

%DAT2TP - Turning points from a signal.

tp = dat2tp(x);    % Calculate turing points
tp5 = dat2tp(x,5); % TP, rainflow filter with threshold 5. Needs FINDRFC!

plot(x(:,1),x(:,2),tp(:,1),tp(:,2),tp5(:,1),tp5(:,2))

pause

%DAT2TP - Turning points from a signal.

[RFM0,u0] = dat2rfm(x);    % Default parameters
[RFM,u] = dat2rfm(x,0.5,[-25 25 50]);
subplot(1,2,1), cmatplot(u0,u0,RFM0,3)
subplot(1,2,2), cmatplot(u,u,RFM,3)

pause

%RFCFILTER - Rainflow filter a signal.
y = rfcfilter(x,0,1); % Filtered signal y is the turning points of x.
y5 = rfcfilter(x,5);  % This removes all rainflow cycles with range less than 5.

plot(x(:,1),x(:,2),y(:,1),y(:,2),y5(:,1),y5(:,2))

pause

%TP2MM - min-max cycles from TP.
[mM,Mm] = tp2mm(tp);
ccplot({mM Mm})

pause

%TP2RFC - Rainflow cycles from TP.
[rfc] = tp2rfc(tp);
[rfc2] = tp2rfc(tp,1);
def=[]; def.res='cs'; def.time=1;
[rfc3] = tp2rfc(tp,def);
def.asymmetric=1;
[rfc4] = tp2rfc(tp,def);
ccplot({rfc rfc2; rfc3 rfc4})

pause

%TP2ARFC - Asymmetric rainflow cycles from TP.
[arfc,arfc0,res] = tp2arfc(tp);
ccplot({arfc arfc0; [] []})
subplot(2,2,3), ccplot(arfc), hold on, plot(arfc0(:,1),arfc0(:,2),'r.'), hold off
subplot(2,2,4), plot(res)

pause

%%%%
% Cycle matrix
%

%CC2CMAT - Cycle matrix from cycle count.
param = [-20 20 41]; u=levels(param);
Frfc = cc2cmat(param,rfc);
Farfc = cc2cmat(param,arfc);
cmatplot(u,u,{Frfc Farfc},13)

pause

% Smoothing, fixed bandwidth
Frfc1 = cc2cmat(param,rfc,[],1);
Farfc1 = cc2cmat(param,arfc,[],1); % Don't work for ARFC
cmatplot(u,u,{Frfc1 Farfc1},13)

pause

% Smoothing, varying bandwidth
Frfc2 = cc2cmat(param,rfc,[],2);
%Farfc2 = cc2cmat(param,arfc,[],2); % Don't work for ARFC
cmatplot(u,u,{Frfc2},13)

pause

%CMAT2NT - Counting distribution from cycle matrix.
F=Frfc;
NT = cmat2nt(F);
 
%NT2CMAT - Cycle matrix from counting distribution.
FF = nt2cmat(NT);
cmatplot({F NT FF},13)
sum(sum(abs(FF-F))) % Shall be zero

pause

%CMAT2RMCMAT - Convert cycle matrix from min-max to range-mean (PJ)
Frm = cmat2rmcmat(F);
cmatplot({F Frm},13)

pause

%RMCMAT2CMAT - Convert cycle matrix from range-mean to min-max (PJ)
F1 = rmcmat2cmat(Frm);
cmatplot({F Frm F1},13)
sum(sum(abs(F1-F))) % Shall be zero

pause

%SMOOTHCMAT - Smooth a cycle matrix. (PJ)
Frfc1 = smoothcmat(Frfc,1); % Smoothing, fixed bandwidth
Frfc2 = smoothcmat(Frfc,2); % Smoothing, varying bandwidth
cmatplot(u,u,{Frfc1 Frfc2},13)

pause

%TP2LC - Level crossings from a cycle count. (WAT/PJ)
lc1 = tp2lc(tp);
subplot(2,2,1), plotlc(lc1)

%CC2LC - Level crossings from a cycle count. (WAT/PJ)
lc2 = cc2lc(rfc);
subplot(2,2,2), plotlc(lc2)

sum(abs(lc1-lc2)) % Shall be zero

%CMAT2LC - Level crossings from cycle matrix. (PJ)
lc3 = cmat2lc(param,Frfc);
subplot(2,2,3), plotlc(lc3)

%NT2LC - Level crossings from counting distribution. (FAT/PJ)
lc4 = nt2lc(param,NT);
subplot(2,2,4), plotlc(lc4)

sum(abs(lc3-lc4)) % Shall be zero

pause

%CC2AMP - Amplitude histogram from a cycle count. (WAT/PJ)
amp = cc2amp(rfc);
[N,X] = hist(amp,((u-u(1))+(u(2)-u(1))/2)/2);
amp_hist1=[X' N'];

%CMAT2AMP - Amplitude histogram from cycle matrix. (PJ)
amp_hist2 = cmat2amp(param,F);

subplot(2,1,1), bar(amp_hist1(:,1),amp_hist1(:,2))
subplot(2,1,2), bar(amp_hist2(:,1)+0.25,amp_hist2(:,2))

lsplot({amp_hist2})
pause

%%%%
% Discrete loads
%

clf

%DAT2DTP - Discrete TP from load. (PJ)
param = [-20 20 41];  u=levels(param);
[dtp1,u1,tp1] = dat2dtp(param,x);
[dtp2,u2,tp2] = dat2dtp(param,x,5);
plot(x(:,1),x(:,2),tp1(:,1),tp1(:,2),tp2(:,1),tp2(:,2)), grid

pause

%CC2DCC - Discretization of a cycle count. (WAT)
dcc = cc2dcc(param,rfc);
Adcc = cc2dcc(param,arfc);

%DTP2RFM - Rainflow matrix from discrete TP. (PJ)
n=param(3);
[RFM,RFM0,res] = dtp2rfm(dtp1(:,2),n);

%DTP2ARFM - Asymmetric rainflow matrix from discrete TP. (PJ)
[ARFM,ARFM0,Ares] = dtp2arfm(dtp1(:,2),n);

%DCC2CMAT - Histogram matrix from discrete class indices. (WAT)
Frfc = dcc2cmat(dcc,n);
Farfc = dcc2cmat(Adcc,n);

% Different methods for discretization, results should be different
cmatplot(u,u,{RFM ARFM; Frfc Farfc},13)

pause



%%%%
% Extrapolation & Smoothing
%

%Create a MCTP model

n = 64;
param = [-1 1 n]; uu=levels(param);
paramD = [1 n n];
[G,Gh] = mktestmat(param,[-0.2 0.2], 0.15,1);
lcG = cmat2lc(paramD,G);

% Simulate a sample path
N = 1000;
xD = mctpsim({G Gh},2*N);

% Calculate RFM
Frfc = dtp2rfm(xD,n);

% Test rfmextrapolate

[G,Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1);
xD = mctpsim({G Gh},2000);
Frfc = dtp2rfm(xD,64,'CS');
Fest = rfmextrapolate(Frfc,[],1);
Grfc = mctp2rfm({G Gh});
cmatplot({Frfc Fest; Grfc G},4)

pause

% Test tpextrapolate

x = load('sea.dat');
tp = dat2tp(x,0.5);
tpe = tpextrapolate(tp,1,[],1);
clf, plot(tp(:,1),tp(:,2),'b',tpe(:,1),tpe(:,2),'r')
[tpe,Pout,I] = tpextrapolate(tp,1,[],2);
clf, plot(tpe(:,1),tpe(:,2),'b',tpe(I.min,1),tpe(I.min,2),'g*',tpe(I.max,1),tpe(I.max,2),'g*')

% Test cmat2extralc

[G,Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1);
xD = mctpsim({G Gh},2000);
Frfc = dtp2rfm(xD,64,'CS');
[lcEst,Est] = cmat2extralc([-1 1 64],Frfc,[-0.4 0.4]);
lcG = cmat2lc([-1 1 64],G/sum(sum(G)));
lcF = cmat2lc([-1 1 64],Frfc);
clf, semilogx(1000*lcG(:,2),lcG(:,1),lcF(:,2),lcF(:,1),lcEst.lc(:,2),lcEst.lc(:,1))

param=[-1 1 64]; u=0.4;
[lcEst,Est] = cmat2extralc(param,Frfc,[-u u]); %,method,plotflag)
lcEst = cmat2extralc(param,Frfc,[-u u],[],0); %,method,plotflag)
[lcEst,Est,R,MSE] = cmat2extralc(param,Frfc,[-u u],'exp,ls1',1); %,method,plotflag)
[lcEst,Est,R,MSE] = cmat2extralc(param,Frfc,[-u u],'exp,wls4',1); %,method,plotflag)
[lcEst,Est,R,MSE] = cmat2extralc(param,Frfc,[-u u],'ray,ml',1); %,method,plotflag)

pause

%Test fitgenpar_mld

data0 = rndgenpar(0.3,1,0,200,1);
x_ML = fzero(@(x)fitgenparml(x,data0),0);
[f,k_ML,s_ML] = fitgenparml(x_ML,data0) % Estimates k_ML and s_ML
data1 = floor(data0*10)/10;
x=(0:0.1:(max(data1)+0.1))';
N = histc(data1+0.05,x);
x_MLD = fzero(@(y)fitgenpar_mld(y,[x N]),0);
[f,k_MLD,s_MLD] = fitgenpar_mld(x_MLD,[x N]) % Estimates k_ML and s_ML

%Test lc2rfmextreme

u = (-5:0.2:5)'; lc = [u exp(-u.^2/2)];
[Frfc,u,Nrfc] = lc2rfmextreme(lc);
cmatplot(u,u,{Frfc Nrfc},3)

pause

%Test cmatcombine

F1 = triu(ones(8),0)
F2 = 2*F1
[F,Lim,FF1,FF2]=cmatcombine(F1,F2,2)
Lim=[]; Lim.range=2; Lim.min=4; Lim.max=4;
[F,Lim,FF1,FF2]=cmatcombine(F1,F2,Lim)

%Test extralc

S = jonswap;
x = spec2sdat(S,100000,0.1,[],'random');
lc = dat2lc(x); s = std(x(:,2));
[lcEst,Est] = extralc(lc,s*[-2 2]);
[lcEst,Est] = extralc(lc,s*[-2 2],'exp,ml');
[lcEst,Est] = extralc(lc,s*[-2 2],'ray,ml');

pause

%Test cmatresamp

F = round(5*triu(rand(4,4),1))
FF = cmatresamp(F)

%%%%%
% End of test

echo('off')
