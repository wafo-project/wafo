%% CHAPTER2 Modelling random loads and stochastic waves
%
% CALL: chapter2
%
% Chapter2 contains the commands used in Chapter 2 of the tutorial and
% present some tools for analysis of random functions with
% respect to their correlation, spectral and distributional properties.
% The presentation is divided into three examples: 
%
% Example1 is devoted to estimation of different parameters in the model.
% Example2 deals with spectral densities and
% Example3 presents the use of WAFO to simulate samples of a Gaussian
% process.
%
% Some of the commands are edited for fast computation. 
% Each set of commands is followed by a 'pause' command.
% Type 'pause off' or ¨pause on' to disable or activa 'pause'.
%
% The figures can be printed in pdf format by setting the parameter
% 'printing' to 1. 

% Tested on Matlab 5.3, 7.10, 8.1, 8.6
% History
% Revised help text and added automatic printing by GL 12-28-2015.
% Revised by Georg Lindgren march 2011 for use with Tutorial 2.5 
% and sept 2009 for WAFO ver 2.5 on Matlab 7.1
% Revised pab sept2005
% Added sections -> easier to evaluate using cellmode evaluation.
% Revised pab Dec2004
% Created by GL July 13, 2000
% from commands used in Chapter 2
%

start=clock;
pstate =  'off';
pause(pstate)
printing=0;

%% Section 2.1 Introduction and preliminary analysis
%% Example 1: Sea data
%% Observed crossings compared to the expected for Gaussian signals
xx = load('sea.dat');
me = mean(xx(:,2))
sa = std(xx(:,2))
xx(:,2) = xx(:,2) - me;
lc = dat2lc(xx);
plotflag = 2;
lcplot(lc,plotflag,0,sa)
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 1');
    print -dpdf ./bilder/C2_1.pdf 
end    
disp('Block = 1'),pause

%% Average number of upcrossings per time unit
% Next we compute the mean frequency as the average number of upcrossings 
% per time unit of the mean level (= 0); this may require interpolation in the 
% crossing intensity curve, as follows.  
T = max(xx(:,1))-min(xx(:,1))
f0 = interp1(lc(:,1),lc(:,2),0)/T  % zero up-crossing frequency 
disp('Block = 2'),pause

%% Turningpoints and irregularity factor
tp = dat2tp(xx);
fm = length(tp)/(2*T)            % frequency of maxima
alfa = f0/fm                     % approx Tm24/Tm02
disp('Block = 3'),pause

%% Visually examine data
% We finish this section with some remarks about the quality
% of the measured data. Especially sea surface measurements can be
% of poor quality. We shall now check the  quality of the dataset {\tt xx}. 
% It is always good practice to visually examine the data 
% before the analysis to get an impression of the quality, 
% non-linearities and narrow-bandedness of the data.
% First we shall plot the data and zoom in on a specific region. 
% A part of sea data is visualized with the following commands
clf
waveplot(xx,tp,'k-','*',1,1)
axis([0 2 -inf inf])   
wafostamp('','(ER)',0)
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 2');
    print -dpdf ./bilder/C2_2.pdf 
end    
disp('Block = 4'),pause

%% Finding possible spurious points
% However, if the amount of data is too large for visual examinations one
% could use the following criteria to find possible spurious points. One
% must be careful using the criteria for extremevalue analysis, because
% it might remove extreme waves that are OK and not spurious.
dt = diff(xx(1:2,1));
dcrit = 5*dt;
ddcrit = 9.81/2*dt*dt;
zcrit = 0;
[inds, indg] = findoutliers(xx,zcrit,dcrit,ddcrit);
disp('Block = 5'),pause

%% Section 2.2 Frequency Modeling of Load Histories
%% Section 2.2.1 Power spectrum, periodogram
% The Periodogram, also called the poser spectrum of the data,  
% separates the energy in the signal over different frequencies

clf
Lmax = 9500;
S = dat2spec(xx,Lmax);
plotspec(S)
axis([0 5 0 0.7])
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 3');
    print -dpdf ./bilder/C2_3.pdf 
end  
disp('Block = 6'),pause

%% Calculate moments  
[mom momtext]= spec2mom(S,4)
[sa sqrt(mom(1))]
disp('Block = 7'),pause

%% Section 2.2.2 Random functions in Spectral Domain - Gaussian processes
%% Smoothing of spectral estimate 
% By decreasing Lmax the spectrum estimate becomes smoother.

% Example 1, contd.
clf
Lmax0 = 200; Lmax1 = 50;
S1 = dat2spec(xx,Lmax0);
S2 = dat2spec(xx,Lmax1);
plotspec(S1,[],'-.')
hold on
plotspec(S2)
hold off
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 4');
    print -dpdf ./bilder/C2_4.pdf 
end    
disp('Block = 8'),pause

%% Estimated autocovariance
% Obviously knowing the spectrum one can compute the covariance
% function. The following matlab code will compute the covariance for the 
% unimodal spectral density S1 and compare it with estimated 
% covariance of the signal xx.
clf
Lmax = 80;
R1 = spec2cov(S1,1);
Rest = dat2cov(xx,Lmax);
covplot(R1,Lmax,[],'.')
hold on
covplot(Rest)
wafostamp('','(ER)')
hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 5');
    print -dpdf ./bilder/C2_5.pdf 
end    
disp('Block = 9'),pause

%%
% We can see in the figure below that the covariance function corresponding 
% to the spectral density S2 significantly differs from the one estimated 
% directly from data. 
% It can be seen in the figure above that the covariance corresponding to S1 
% agrees much better with the estimated covariance function.
R2 = spec2cov(S2,1);
covplot(R2,Lmax,[],'.')
hold on
covplot(Rest)
wafostamp('','(ER)')
hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 6');
    print -dpdf ./bilder/C2_6.pdf 
end    
disp('Block = 11'),pause

%% Section 2.2.4 Transformed Gaussian models
% We begin with computing skewness and kurtosis
% for the data set xx and compare it with the second order wave approximation
% proposed by Winterstein:
rho3 = skew(xx(:,2))
rho4 = kurt(xx(:,2))
[sk, ku]=spec2skew(S1)
disp('Block = 12'),pause

%% Comparisons of 3 transformations
clf
gh = hermitetr([],[sa sk ku me]);
g  = gh; g(:,2)=g(:,1)/sa;  % Linear transformation
trplot(g)  
[glc, test0 cmax irr gemp] = dat2tr(xx,[],'plotflag',1);
hold on
plot(glc(:,1),glc(:,2),'b-') % Transf. estimated from level-crossings
plot(gh(:,1),gh(:,2),'b-.')  % Hermite Transf. estimated from moments 
hold off
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 7');
    print -dpdf ./bilder/C2_7.pdf 
end    
disp('Block = 13'),pause

%%  Test Gaussianity of a stochastic process.
% TESTGAUSSIAN simulates  e(g(u)-u) = int (g(u)-u)^2 du  for Gaussian processes 
% given the spectral density, S. The result is plotted if test0 is given.
% This is useful for testing if the process X(t) is Gaussian.
% If 95% of TEST1 is less than TEST0 then X(t) is not Gaussian at a 5% level.
% 
% As we see from the figure below: none of the simulated values of test1 is 
% above 1.00. Thus the data significantly departs from a Gaussian distribution. 

N = length(xx);
clf
test1 = testgaussian(S1,[N,50],test0);
wafostamp('','(CR)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 8');
    print -dpdf ./bilder/C2_8.pdf 
end    
disp('Block = 14'),pause

%% Normalplot of data xx
% indicates that the underlying distribution has a "heavy" upper tail and a
% "light" lower tail. 
clf
plotnorm(xx(:,2))
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 9');
    print -dpdf ./bilder/C2_9.pdf 
end    
disp('Block = 15'),pause

%% Section 2.2.5 Spectral densities of sea data
%% Example 2: Different forms of spectra
clf
Hm0 = 7; Tp = 11;
spec = jonswap([],[Hm0 Tp]);
spec.note
disp('Block = 16'),pause

%% Directional spectrum and Encountered directional spectrum
%% Directional spectrum
clf
D = spreading(101,'cos2s',0,[],spec.w,1)
Sd = mkdspec(spec,D)
disp('Block = 17'),pause

%% Encountered directional spectrum
clf
Se = spec2spec(Sd,'encdir',0,10);
plotspec(Se), hold on
plotspec(Sd,1,'--'), hold off
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 10');
    print -dpdf ./bilder/C2_10.pdf 
end    
disp('Block = 18'),pause

%% Frequency spectra
clf
Sd1 =spec2spec(Sd,'freq');
Sd2 = spec2spec(Se,'enc');
plotspec(spec), hold on
plotspec(Sd1,1,'.'),
plotspec(Sd2),
wafostamp('','(ER)')
hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 11');
    print -dpdf ./bilder/C2_11.pdf 
end    
disp('Block = 19'),pause

%% Wave number spectrum
clf
Sk = spec2spec(spec,'k1d')
Skd = spec2spec(Sd,'k1d')
plotspec(Sk), hold on
plotspec(Skd,1,'--'), hold off
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 12');
    print -dpdf ./bilder/C2_12.pdf 
end    
disp('Block = 20'),pause

%% Effect of waterdepth on spectrum
clf
plotspec(spec,1,'--'), hold on
S20 = spec;
S20.S = S20.S.*phi1(S20.w,20);
S20.h = 20;
plotspec(S20),  hold off
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 13');
    print -dpdf ./bilder/C2_13.pdf 
end    
disp('Block = 21'),pause

%% Section 2.3 Simulation of transformed Gaussian process
%% Example 3: Simulation of random sea    
% The reconstruct function replaces the spurious points of seasurface by
% simulated data on the basis of the remaining data and a transformed Gaussian
% process. As noted previously one must be careful using the criteria 
% for finding spurious points when reconstructing a dataset, because
% these criteria might remove the highest and steepest waves as we can see
% in this plot where the spurious points is indicated with a '+' sign:
%
clf
[y, grec] = reconstruct(xx,inds); pause(pstate)
waveplot(y,'-',xx(inds,:),'+',1,1)
axis([0 inf -inf inf])
wafostamp('','(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 14');
    print -dpdf ./bilder/C2_14.pdf 
end    
disp('Block = 22'), pause

% Compare transformation (grec) from reconstructed (y) 
% with original (glc) from (xx)
clf
trplot(g), hold on
plot(gemp(:,1),gemp(:,2))
plot(glc(:,1),glc(:,2),'-.')
plot(grec(:,1),grec(:,2)), hold off 
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 15');
    print -dpdf ./bilder/C2_15.pdf 
end    
disp('Block = 23'),pause

%%
clf
L = 200;
x = dat2gaus(y,grec);
Sx = dat2spec(x,L);
disp('Block = 24'),pause
      
%%
clf
dt = spec2dt(Sx)
Ny = fix(2*60/dt) % = 2 minutes
Sx.tr = grec;
ysim = spec2sdat(Sx,Ny);
waveplot(ysim,'-')
wafostamp('','(CR)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 16');
    print -dpdf ./bilder/C2_16.pdf 
end    
disp('Block = 25'),pause
 
%% Estimated spectrum compared to Torsethaugen spectrum
clf
Tp = 1.1;
H0 = 4*sqrt(spec2mom(S1,1))
St = torsethaugen([0:0.01:5],[H0  2*pi/Tp]);
plotspec(S1)
hold on
plotspec(St,'-.')
axis([0 6 0 0.4])
wafostamp('','(ER)'); hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 17');
    print -dpdf ./bilder/C2_17.pdf 
end    
disp('Block = 26'),pause

%%
clf
Snorm = St;
Snorm.S = Snorm.S/sa^2;
dt = spec2dt(Snorm)
disp('Block = 27'),pause

%%
clf
[Sk Su] = spec2skew(St);
sa = sqrt(spec2mom(St,1));
gh = hermitetr([],[sa sk ku me]);
Snorm.tr = gh;
disp('Block = 28'),pause

%% Transformed Gaussian model compared to Gaussian model
clf
dt = 0.5;
ysim_t = spec2sdat(Snorm,240,dt);
xsim_t = dat2gaus(ysim_t,Snorm.tr);
disp('Block = 29'),pause

%% Compare
% In order to compare the Gaussian and non-Gaussian models we need to scale  
% \verb+xsim_t+ %{\tt xsim$_t$} 
% to have the same first spectral moment as 
% \verb+ysim_t+,  %{\tt ysim$_t$},  Since the process xsim_t has variance one
% which will be done by the following commands. 
clf
xsim_t(:,2) = sa*xsim_t(:,2);
waveplot(xsim_t,ysim_t,5,1,sa,4.5,'r.','b')
wafostamp('','(CR)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C2 18');
    print -dpdf ./bilder/C2_18.pdf 
end    
disp('Block = 30, Last block'),pause

disp('Elapsed time')
etime(clock,start)
