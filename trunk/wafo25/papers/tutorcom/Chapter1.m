%% CHAPTER1 demonstrates some applications of WAFO
%
% CHAPTER1 gives an overview through examples some of the capabilities of
% WAFO. WAFO is a toolbox of Matlab routines for statistical analysis and
% simulation of random waves and loads.
%
% The commands are edited for fast computation.
% Each set of commands is followed by a 'pause' command.
% Type 'pause off' to disable them.

% Tested on Matlab 5.3, 7.0
% History
% Revised pab sept2005
%  Added sections -> easier to evaluate using cellmode evaluation.
% Revised pab Dec 2004
% Added support for publish.m command in matlab R14
% Created by GL July 13, 2000
% from commands used in Chapter 1 of the tutorial

pstate = 'off'

%% Section 1.4 Some applications of WAFO

%% Section 1.4.1 Simulation from spectrum, estimation of spectrum 
%% Simulation of the sea surface from spectrum
%The following code generates 200 seconds of data sampled with 10Hz from
%the Torsethaugen spectrum
Hm0 = 6;
Tp  = 8;
S1=torsethaugen([],[Hm0 Tp],1);
clf
dt = 0.1;
N = 2000;
xs=spec2sdat(S1,N,dt);

clf
waveplot(xs,'-')
wafostamp([],'(ER)')
pause(pstate)

%% Estimation of spectrum 
%A common situation is that one wants to estimate the spectrum for wave
%measurements. The following code simulate 20 minutes signal sampled at 4Hz
%and compare the spectral estimate with the original Torsethaugen spectum.
clf
xs=spec2sdat(S1,[20*60*4 1],0.25);
Sest = dat2spec(xs,400);
plotspec(Sest,1,'--'), hold on
plotspec(S1,1), hold off
axis([0 3 0 5])
wafostamp([],'(ER)')
pause(pstate)


%% Section 1.4.2 Probability distributions of wave characteristics.
%% Probability distribution of wave trough period
%WAFO gives the possibility of computing the exact probability
%distributions for a number of characteristics given a spectral density.
%In the following example we study the trough period extracted from the
%time series and compared with the theoretical density computed with exact
%spectrum, S1, and the estimated spectrum, Sest.

clf
[T, index] = dat2wa(xs,0,'d2u');
histgrm(T,25,1,1), hold on
dtyex = spec2tpdf(S1,[],'Tt',[0 10 51],0,3);
dtyest = spec2tpdf(Sest,[],'Tt',[0 10 51],0,3);
pdfplot(dtyex)
pdfplot(dtyest,'-.')
axis([0 10 0 0.35]), hold off
wafostamp([],'(ER)')
pause(pstate)

%% Section 1.4.3 Directional spectra
%Here are a few lines of code, which produce directional spectra 
%with frequency independent and frequency dependent spreading.
clf
D1 = spreading(101,'cos',pi/2,15,[],0); % frequency independent
D12 = spreading(101,'cos',0,15,S1.w,1); % frequency dependent
SD1 = mkdspec(S1,D1);
SD12 = mkdspec(S1,D12);
plotspec(SD1,1), hold on, plotspec(SD12,1,'-.'); hold off
wafostamp([],'(ER)')
pause(pstate)


%% 3D Simulation of the sea surface 
% The simulations show that frequency dependent spreading leads to
% much more irregular surface so the orientation of waves is less
% transparent compared to the frequency independent case.

% Frequency independent spreading
Y1=seasim(SD1,2^8,2^8,1,0.5,0.5,0.25,2,1);
wafostamp([],'(ER)')
pause(pstate)
%%
% Frequency dependent spreading
Y12=seasim(SD12,2^8,2^8,1,0.5,0.5,0.25,2,1);
wafostamp([],'(ER)')
pause(pstate)

%% Section 1.4.4 Fatigue, Load cycles and Markov models.
%% Switching Markow chain of turningpoints 
% Here the Markov approximation for computing the intensity of
% rainflowcycles for the Gaussian model with spectrum S1
clf
frfc=spec2cmat(S1,[],'rfc',[],[-6 6 61]);
pdfplot(frfc);
hold on
tp=dat2tp(xs);
rfc=tp2rfc(tp);
plot(rfc(:,2),rfc(:,1),'.')
wafostamp([],'(ER)')
hold off

%% 1.4.5 Statistical extreme value analysis
% The WAFO-toolbox contains almost 300 routines for general statistical analysis, 
% description, plotting, and simulation. In Chapter 5 we describe some routines which 
% are particularly important for wave and fatigue analysis, related to statistics of extremes. These are
% based on the generalized extreme value (GEV) and generalized Pareto distribution (GPD),
% combined with the peaks over threshold (POT) method.
% As an example we show an analysis of wave elevation data from the Poseidon plat-
% form in the Japan Sea. Data from about 23 hours of registration are stored in the data set
% yura87, taken with a 1 Hz sampling rate. We ?rst load and plot, in Figure 1.7, part of
% the data and calculate the maximum over 5 minute periods.
clf()
xn=load('yura87.dat'); subplot(211);
plot(xn(1:30:end,1)/3600,xn(1:30:end,2),'.')
title('Water level'), ylabel('m')
yura=xn(1:85500,2);
yura=reshape(yura,300,285);
maxyura=max(yura);
subplot(212)
plot(xn(300:300:85500,1)/3600,maxyura,'.')
xlabel('Time (h)'), ylabel('m')
title('Maximum 5 min water level')

%%
% It is clear from the ?gures that there is a trend in the data, with decreasing spreading
% with time. In Chapter 5 we will deal with that problem; here we make a crude extreme
% value analysis, by ?tting a GEV distribution to the sequence of 5 minute maxima, simply
% by issuing the commands
figure()
phat=fitgev(maxyura,'plotflag',1);

