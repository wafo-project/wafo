%% CHAPTER6 contains the commands used in Chapter 6 of the tutorial
%
% CALL:  Chapter6
% 
% Some of the commands are edited for fast computation. 
% Each set of commands is followed by a 'pause' command.
% 

% Tested on Matlab 5.3, 7.10
% History
% Revised by Georg Lindgren march 2011 for use with Tutorial 2.5 and 
% sept 2009 for WAFO ver 2.5 on Matlab 7.1
% Added Return values by GL August 2008
% Revised pab sept2005
% Added sections -> easier to evaluate using cellmode evaluation.
% Created by GL July 13, 2000
% from commands used in Chapter 5
%

%% Chapter 6 Extreme value analysis
start=clock;
pstate = 'off';
pause(pstate)

%% Section 6.1 Weibull and Gumbel papers
% Significant wave-height data on Weibull paper,
clf
Hs = load('atlantic.dat');
wei = plotweib(Hs)
wafostamp([],'(ER)')
disp('Block = 1'),pause

%%
% Significant wave-height data on Gumbel paper,
clf
gum=plotgumb(Hs)
wafostamp([],'(ER)')
  disp('Block = 2'),pause

%%
% Significant wave-height data on Normal probability paper,
plotnorm(log(Hs),1,0);
wafostamp([],'(ER)')
disp('Block = 3'),pause

%%
% Return values in the Gumbel distribution
clf
T=1:100000;
sT=gum(2) - gum(1)*log(-log(1-1./T));
semilogx(T,sT), hold
N=1:length(Hs); Nmax=max(N);
plot(Nmax./N,sort(Hs,'descend'),'.')
title('Return values in the Gumbel model')
xlabel('Return period')
ylabel('Return value') 
wafostamp([],'(ER)')
disp('Block = 4'),pause

%% Section 6.2 Generalized Pareto and Extreme Value distributions
%% Section 6.2.1 Generalized Extreme Value distribution

% Empirical distribution of significant wave-height with estimated 
% Generalized Extreme Value distribution,
gev=fitgev(Hs,'plotflag',2)
wafostamp([],'(ER)')
disp('Block = 5a'),pause

clf
x = linspace(0,14,200);
plotkde(Hs,[x;pdfgev(x,gev)]')
disp('Block = 5b'),pause

% Analysis of yura87 wave data. 
% Wave data interpolated (spline) and organized in 5-minute intervals
% Normalized to mean 0 and std = 1 to get stationary conditions. 
% maximum level over each 5-minute interval analysed by GEV
xn  = load('yura87.dat');
XI  = 0:0.25:length(xn);
N   = length(XI); N = N-mod(N,4*60*5); 
YI  = interp1(xn(:,1),xn(:,2),XI(1:N),'spline');
YI  = reshape(YI,4*60*5,N/(4*60*5)); % Each column holds 5 minutes of 
                                     % interpolated data.
Y5  = (YI-ones(1200,1)*mean(YI))./(ones(1200,1)*std(YI)); 
Y5M = max(Y5);
Y5gev = fitgev(Y5M,'method','mps','plotflag',2)
wafostamp([],'(ER)')
disp('Block = 6'),pause

%% Section 6.2.2 Generalized Pareto distribution

% Exceedances of significant wave-height data over level 3,
gpd3 = fitgenpar(Hs(Hs>3)-3,'plotflag',1);
wafostamp([],'(ER)')

%%
figure
% Exceedances of significant wave-height data over level 7,
gpd7 = fitgenpar(Hs(Hs>7),'fixpar',[nan,nan,7],'plotflag',1);
wafostamp([],'(ER)')
disp('Block = 6'),pause

%% 
%Simulates 100 values from the GEV distribution with parameters (0.3, 1, 2), then estimates the
%parameters using two different methods and plots the estimated distribution functions together
%with the empirical distribution.
Rgev = rndgev(0.3,1,2,1,100);
gp = fitgev(Rgev,'method','pwm');
gm = fitgev(Rgev,'method','ml','start',gp.params,'plotflag',0);

x=sort(Rgev);
plotedf(Rgev,gp,{'-','r-'});
hold on
plot(x,cdfgev(x,gm),'--')
hold off
wafostamp([],'(ER)')
disp('Block =7'),pause

%%
% Similarly for the GPD distribution;
Rgpd = rndgenpar(0.4,1,0,1,100);
plotedf(Rgpd);
hold on
gp = fitgenpar(Rgpd,'method','pkd','plotflag',0);
x=sort(Rgpd);
plot(x,cdfgenpar(x,gp))
% gm = fitgenpar(Rgpd,'method','mom','plotflag',0);
% plot(x,cdfgenpar(x,gm),'g--')
gw = fitgenpar(Rgpd,'method','pwm','plotflag',0);
plot(x,cdfgenpar(x,gw),'g:')
gml = fitgenpar(Rgpd,'method','ml','plotflag',0);
plot(x,cdfgenpar(x,gml),'--')
gmps = fitgenpar(Rgpd,'method','mps','plotflag',0);
plot(x,cdfgenpar(x,gmps),'r-.')
hold off
wafostamp([],'(ER)')
disp('Block = 8'),pause

%%
% Return values for the GEV distribution
T = logspace(1,5,10);
[sT, sTlo, sTup] = invgev(1./T,Y5gev,'lowertail',false,'proflog',true);

%T = 2:100000;
%k=Y5gev.params(1); mu=Y5gev.params(3); sigma=Y5gev.params(2);
%sT1 = invgev(1./T,Y5gev,'lowertail',false);
%sT=mu + sigma/k*(1-(-log(1-1./T)).^k);
clf
semilogx(T,sT,T,sTlo,'r',T,sTup,'r'), hold
N=1:length(Y5M); Nmax=max(N);
plot(Nmax./N,sort(Y5M,'descend'),'.')
title('Return values in the GEV model')
xlabel('Return priod')
ylabel('Return value') 
grid on 
disp('Block = 9'),pause

%% Section 6.3 POT-analysis

% Estimated expected exceedance over level u as function of u.
clf
plotreslife(Hs,'umin',2,'umax',10,'Nu',200);
wafostamp([],'(ER)')
disp('Block = 10'),pause

%%
% Estimated distribution functions of monthly maxima 
%with the POT method (solid),
% fitting a GEV (dashed) and the empirical distribution.

% POT- method
gpd7 = fitgenpar(Hs(Hs>7)-7,'method','pwm','plotflag',0);
khat = gpd7.params(1);
sigmahat = gpd7.params(2);
muhat = length(Hs(Hs>7))/(7*3*2);
bhat = sigmahat/muhat^khat;
ahat = 7-(bhat-sigmahat)/khat;
x = linspace(5,15,200);
plot(x,cdfgev(x,khat,bhat,ahat))
disp('Block = 11'),pause

%%
% Since we have data to compute the monthly maxima over 42 months 
% we can also try to fit a GEV distribution directly:
mm = zeros(1,41);
for i=1:41
  mm(i)=max(Hs(((i-1)*14+1):i*14)); 
end

gev=fitgev(mm);
hold on
plotedf(mm)
plot(x,cdfgev(x,gev),'--')
hold off
wafostamp([],'(ER)')
disp('Block = 12'),pause

%%
% Section 6.3.3 Declustering

Nmin=7; % minimum number of extremes
Tmin=5; % minimum distance between extremes
Tb = 15; % block period
xn = load('sea.dat');
timeSpan = (xn(end,1)-xn(1,1))/60; % in minutes
dt = xn(2,1)-xn(1,1); % in seconds
tc = dat2tc(xn);
umin = median(tc(tc(:,2)>0,2));
Ie0 = findpot(tc, 0.9*umin, Tmin);
Ev = sort(tc(Ie0,2));
Ne = numel(Ev)
if Ne>Nmin && Ev(Ne-Nmin)>umin,
   umax = Ev(Ne-Nmin);
else
  umax = umin;
end

Nu = floor((umax-umin)/0.025)+1;
u = linspace(umin,umax,Nu);
mrl = reslife(Ev, 'u',u);
umin0 = umin;
for io= numel(mrl.data):-1:1,
    CI = mrl.dataCI(io:end,:);
    if ~(max(CI(:,1))<= mrl.data(io) & mrl.data(io)<=min(CI(:,2))),
        umin0 = mrl.args(io);
    break;
    end
end
[di, threshold, ok_u] = disprsnidx(tc(Ie0,:), 'Tb', Tb, 'alpha',0.05, 'u',u);
disp('Block = 13'), pause

%% 
% Plotting

figure(1); plot(di)
vline(threshold) % Threshold from dispersion index
vline(umin0,'g') % Threshold from mean residual life plot
figure(2); plot(mrl)
vline(threshold) % Threshold from dispersion index
vline(umin0,'g') % Threshold from mean residual life plot

% threshold around 1.2 seems appropriate

Ie = findpot(tc, threshold, Tmin);
lambda = numel(Ie)/timeSpan; % # Y> threshold per minute
varLambda = lambda*(1-(dt/60)*lambda)/timeSpan;
stdLambd = sqrt(varLambda)
Ev = tc(Ie,2);
phat = fitgenpar(Ev, 'fixpar',[nan,nan,threshold], 'method','mps');

figure(3); phat.plotfitsumry() % check fit to data
disp('Block = 14'), pause

%%
% Calculate return period and return level
%Tr = 60; %Return period in minutes
%[xr,xrlo,xrup] = invgenpar(1./(lambda*Tr),phat,'lowertail',false,...
%    'alpha', 0.05); % 60 minutes return level + 95%CI
%[xr,xrlo5,xrup5] = invgenpar(1./(lambda*Tr),phat,'lowertail',false,...
%    'alpha', 0.5); % 60 minutes return level + 50%CI

Tr = 3*60; %Return period in minutes
[xr,xrlo,xrup] = invgenpar(1./(lambda*Tr),phat,'lowertail',false,...
    'alpha', 0.05); % 180 minutes return level + 95%CI
[xr,xrlo5,xrup5] = invgenpar(1./(lambda*Tr),phat,'lowertail',false,...
    'alpha', 0.5); % 180 minutes return level + 50%CI

disp('Block 15, Final block')

disp('Elapsed time')
etime(clock,start)
