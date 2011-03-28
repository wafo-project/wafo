%% CHAPTER1 demonstrates some applications of WAFO
%
% CHAPTER1 gives an overview through examples some of the capabilities of
% WAFO. WAFO is a toolbox of Matlab routines for statistical analysis and
% simulation of random waves and loads.
%
% The commands are edited for fast computation.
% Each set of commands is followed by a 'pause' command.
% Type 'pause off' to disable them.

% Tested on Matlab 5.3, 7.10
% History
% Revised by Georg Lindgren March 2011 for use with Tutorial 2.5 and 
% sept 2009 for WAFO ver 2.5 on Matlab 7.1
% Revised pab sept2005
% Added sections -> easier to evaluate using cellmode evaluation.
% Revised pab Dec 2004
% Added support for publish.m command in matlab R14
% Created by GL July 13, 2000
% from commands used in Chapter 1 of the tutorial

start=clock;
pstate = 'off'
pause(pstate)

%% Section 1.4 Some applications of WAFO

%% Section 1.4.1 Simulation from spectrum, estimation of spectrum 
%% Simulation of the sea surface from spectrum
% The following code generates 200 seconds of data sampled with 10Hz from
% the Torsethaugen spectrum
Hm0 = 6;
Tp  = 8;
plotflag = 1;
S1=torsethaugen([],[Hm0 Tp],plotflag);
disp('Block = 1'),pause

%%
dt = 0.1;
N = 2000;
xs=spec2sdat(S1,N,dt);

clf
waveplot(xs,'-')
%wafostamp('','(ER)')
disp('Block = 2'),pause

%% Estimation of spectrum 
% A common situation is that one wants to estimate the spectrum for wave
% measurements. The following code simulate 20 minutes signal sampled at 4Hz
% and compare the spectral estimate with the original Torsethaugen spectum.
clf
plotflag = 1;
Fs = 4; 
dt = 1/Fs;
N  = fix(20*60*Fs);
xs   = spec2sdat(S1,N,dt);
Sest = dat2spec(xs,400)
plotspec(S1,plotflag), hold on
plotspec(Sest,plotflag,'--'), hold off
axis([0 3 0 5]) % This may depend on the simulation
%wafostamp('','(ER)')
disp('Block = 3'),pause

%% Section 1.4.2 Probability distributions of wave characteristics.
%% Probability distribution of wave trough period
% WAFO gives the possibility of computing the exact probability
% distributions for a number of characteristics given a spectral density.
% In the following example we study the trough period extracted from the
% time series and compared with the theoretical density computed with exact
% spectrum, S1, and the estimated spectrum, Sest.

clf
NIT = 3
paramt = [0 10 51];
dtyex = spec2tpdf(S1,[],'Tt',paramt,0,NIT);
dtyest = spec2tpdf(Sest,[],'Tt',paramt,0,NIT);
[T, index] = dat2wa(xs,0,'d2u');
histgrm(T,25,1,1), hold on
pdfplot(dtyex)
pdfplot(dtyest,'-.')
axis([0 10 0 0.35]), hold off
%wafostamp('','(ER)')
disp('Block = 4'),pause

%% Section 1.4.3 Directional spectra
% Here are a few lines of code, which produce directional spectra 
% with frequency independent and frequency dependent spreading.
clf
plotflag = 1
Nt = 101;   % number of angles
th0 = pi/2; % primary direction of waves
Sp  = 15;   % spreading parameter
D1 = spreading(Nt,'cos',th0,Sp,[],0); % frequency independent
D12 = spreading(Nt,'cos',0,Sp,S1.w,1); % frequency dependent
SD1 = mkdspec(S1,D1);
SD12 = mkdspec(S1,D12);
plotspec(SD1,plotflag), hold on, plotspec(SD12,plotflag,'-.'); hold off
wafostamp('','(ER)')
disp('Block = 5'),pause

%% 3D Simulation of the sea surface 
% The simulations show that frequency dependent spreading leads to
% much more irregular surface so the orientation of waves is less
% transparent compared to the frequency independent case.
%
% Frequency independent spreading
plotflag = 1; iseed = 1;
Nx = 2^8;Ny = Nx;Nt = 1;dx = 0.5; dy = dx; dt = 0.25; fftdim = 2;
randn('state',iseed)
Y1 = seasim(SD1,Nx,Ny,Nt,dx,dy,dt,fftdim,plotflag);
wafostamp('','(ER)')
axis('fill')
disp('Block = 6'),pause

%%
% Frequency dependent spreading
randn('state',iseed)
Y12 = seasim(SD12,Nx,Ny,Nt,dx,dy,dt,fftdim,plotflag);
wafostamp('','(ER)')
axis('fill')
disp('Block = 7'),pause

%% Estimation of directional spectrum
%  The figure is not shown in the Tutorial

 Nx = 3; Ny = 2; Nt = 2^12; dx = 10; dy = 10;dt = 0.5;
 F  = seasim(SD12,Nx,Ny,Nt,dx,dy,dt,1,0);  
 Z  = permute(F.Z,[3 1 2]);
 [X,Y] = meshgrid(F.x,F.y);
 N = Nx*Ny;
 types = repmat(sensortypeid('n'),N,1);
 bfs   = ones(N,1);
 pos   = [X(:),Y(:),zeros(N,1)];
 h = inf;
 nfft = 128;
 nt = 101;
 SDe = dat2dspec([F.t Z(:,:)],[pos types,bfs],h,nfft,nt);
plotspec(SDe), hold on
plotspec(SD12,'--'), hold off
disp('Block = 8'),pause

%% Section 1.4.4 Fatigue, Load cycles and Markov models.
%% Switching Markow chain of turningpoints 
% In fatigue applications the exact sample path is not important, but
% only the tops and bottoms of the load, called the sequence of turning
% points (TP). From the turning points one can extract load cycles, from
% which damage calculations and fatigue life predictions can be
% performed.
%
% The matlab commands below computes the intensity of rainflowcycles for 
% the Gaussian model with spectrum S1 using the Markow approximation. 
% The rainflow cycles found in the simulated load signal are shown in the 
% figure.
clf
paramu = [-6 6 61];
frfc=spec2cmat(S1,[],'rfc',[],paramu);
pdfplot(frfc);
hold on
tp = dat2tp(xs);
rfc = tp2rfc(tp);
plot(rfc(:,2),rfc(:,1),'.')
wafostamp('','(ER)')
hold off
disp('Block = 9'),pause

%% Section 1.4.5 Extreme value statistics
% Plot of yura87 data
clf
xn=load('yura87.dat'); subplot(211); 
plot(xn(1:30:end,1)/3600,xn(1:30:end,2),'.')
title('Water level','FontSize',16)
ylabel('(m)','FontSize',16)
set(gca,'FontSize',14)
% Formation of 5 min maxima
yura=xn(1:85500,2);
yura=reshape(yura,300,285);
maxyura=max(yura);
subplot(212)
plot(xn(300:300:85500,1)/3600,maxyura,'.')
xlabel('Time (h)','FontSize',16)
ylabel('(m)','FontSize',16)
title('Maximum 5 min water level','FontSize',16)
set(gca,'FontSize',14)
disp('Block = 10'),pause

%% Estimation of GEV for yuramax
clf
phat=fitgev(maxyura,'plotflag',1);
disp('Block = 11, Last block')
disp('Elapsed time')
etime(clock,start)

