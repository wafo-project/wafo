%% CHAPTER3  Demonstrates empirical distributions of wave characteristics
%
% CALL: chapter3
%
% Chapter3 contains the commands used in Chapter3 in the tutorial.
% 
% Some of the commands are edited for fast computation. 
% Each set of commands is followed by a 'pause' command.
%
% The figures can be printed in pdf format by setting the parameter
% 'printing' to 1. %

% Tested on Matlab 5.3, 7.10, 8.1, 8.6
% History
% Revised by Georg Lindgren march 2011 for use with Tutorial 2.5 and 
% sept 2009 for WAFO ver 2.5 on Matlab 7.1
% Revised pab sept2005
% Added sections -> easier to evaluate using cellmode evaluation.
% Revised by pab Feb 2005
% -updated calls to kdetools+spec2XXpdf programs
% Created by GL July 12, 2000
% from commands used in Chapter 3, written by IR
%

start=clock;
pstate = 'off';
speed = 'fast';
%speed = 'slow'
pause(pstate)
printing=0;

%% Section 3.2 Estimation of wave characteristics from data
%% Section 3.2.1 Wave period 
%Example 1

xx = load('sea.dat');
xx(:,2) = detrend(xx(:,2));
rate = 8;
Tcrcr = dat2wa(xx,0,'c2c','tw',rate);
Tc = dat2wa(xx,0,'u2d','tw',rate);
disp('Block = 1'), pause

%% Histogram of crestperiod compared to the kernel density estimate
clf
mean(Tc)
max(Tc)
t = linspace(0.01,8,200);
kopt = kdeoptset('L2',0);
tic
ftc1 = kde(Tc,kopt,t);
toc
pdfplot(ftc1); hold on
histgrm(Tc,[],[],1)
axis([0 8 0 0.5])
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 1');
    print -dpdf ./bilder/C3_1.pdf 
end    
disp('Block = 2'), pause
%%
clf
tic
kopt.inc = 128;
ftc2 = kdebin(Tc,kopt);
toc
pdfplot(ftc2,'-.');
title('Kernel Density Estimates')
hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 2');
    print -dpdf ./bilder/C3_2.pdf 
end    
disp('Block = 3'), pause

%% Section 3.2.2 Extreme waves - model check: the highest and steepest wave
%clf
clf
method = 0;
rate = 8;
[S, H, Ac, At, Tcf, Tcb, z_ind, yn] = ...
   dat2steep(xx,rate,method);
disp('Block = 4'), pause

%clf
[Smax indS]=max(S)
[Amax indA]=max(Ac)
spwaveplot(yn,[indA indS],'k.')
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 3');
    print -dpdf ./bilder/C3_3.pdf 
end    
disp('Block = 5'), pause

%% Does the highest wave contradict a transformed Gaussian model?
clf
inds1 = (5965:5974)'; % points to remove
Nsim = 10;
[y1, grec1, g2, test, tobs, mu1o, mu1oStd] = ...
   reconstruct(xx,inds1,Nsim); pause(pstate)
spwaveplot(y1,indA-10)
hold on
plot(xx(inds1,1),xx(inds1,2),'+')
lamb = 2.;
muLstd = tranproc(mu1o-lamb*mu1oStd,fliplr(grec1));
muUstd = tranproc(mu1o+lamb*mu1oStd,fliplr(grec1));
plot (y1(inds1,1), [muLstd muUstd],'b-')
axis([1482 1498 -1 3]), 
wafostamp([],'(ER)'); hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 4');
    print -dpdf ./bilder/C3_4.pdf 
end    
disp('Block = 6'), 
pause

%% Expected value (solid) compared to data removed
clf
plot(xx(inds1,1),xx(inds1,2),'+'), hold on
mu = tranproc(mu1o,fliplr(grec1));
plot(y1(inds1,1), mu), hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 5');
    print -dpdf ./bilder/C3_5.pdf 
end    
disp('Block = 7'), pause

%% Section 3.2.3 Crest height PDF
% Transform data so that kde works better
clf
L2 = 0.6;
plotnorm(Ac.^L2); hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 6');
    print -dpdf ./bilder/C3_6.pdf 
end    
pause
%%
%
clf
fac = kde(Ac,{'L2',L2},linspace(0.01,3,200));
pdfplot(fac);
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 7');
    print -dpdf ./bilder/C3_7.pdf 
end    
simpson(fac.x{1},fac.f)
disp('Block = 8'), pause

%% Empirical crest height CDF
clf
Fac = flipud(cumtrapz(fac.x{1},flipud(fac.f)));
Fac = [fac.x{1} 1-Fac];
Femp = plotedf(Ac,Fac);
axis([0 2 0 1])
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 8');
    print -dpdf ./bilder/C3_8.pdf 
end    
disp('Block = 9'), pause

%% Empirical crest height CDF compared to a Transformed Rayleigh approximation
clf
facr = trraylpdf(fac.x{1},'Ac',grec1);
Facr = cumtrapz(facr.x{1},facr.f);
hold on
plot(facr.x{1},Facr,'.')
axis([1.25 2.25 0.95 1])
wafostamp([],'(ER)'); hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 9');
    print -dpdf ./bilder/C3_9.pdf 
end    
disp('Block = 10'), pause

%% Section 3.2.4 Joint pdf of crest period and crest height
clf
kopt2 = kdeoptset('L2',0.5,'inc',256);
Tc = Tcf+Tcb;
fTcAc = kdebin([Tc Ac],kopt2);
fTcAc.labx={'Tc [s]'  'Ac [m]'}
pdfplot(fTcAc);
hold on
plot(Tc,Ac,'k.')
hold off
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 10');
    print -dpdf ./bilder/C3_10.pdf 
end    
disp('Block = 11'), pause

%% Section 3.3 Explicit results - parametric models 
% We investigare three common approximations to wave period and 
% amplitude distributions
%% Section 3.3.1 The average wave
% Example 4: Simple wave characteristics obtained from Jonswap spectrum 
clf
S = jonswap([],[5 10]);
[m,  mt]= spec2mom(S,4,[],0);
disp('Block = 12'), pause

clf
spec2bw(S)
[ch Sa2] = spec2char(S,[1 3]);
disp('Block = 13'), pause

%% Section 3.3.2 Explicit approximations of wave distributions 
%% Longuett-Higgins model for Tc and Ac
%clf
t = linspace(0,15,100);
h = linspace(0,6,100);
flh = lh83pdf(t,h,[m(1),m(2),m(3)]);
pdfplot(flh);
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 11');
    print -dpdf ./bilder/C3_11.pdf 
end    
disp('Block = 14'), pause

%% Transformed Longuett-Higgins model for Tc and Ac
clf
[sk, ku ]=spec2skew(S);
sa = sqrt(m(1));
gh = hermitetr([],[sa sk ku 0]);
flhg = lh83pdf(t,h,[m(1),m(2),m(3)],gh);
pdfplot(flhg);
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 12');
    print -dpdf ./bilder/C3_12.pdf 
end    
disp('Block = 15'), pause

%% Cavanie model for Tc and Ac
clf
t = linspace(0,10,100);
h = linspace(0,7,100);
fcav = cav76pdf(t,h,[m(1) m(2) m(3) m(5)],[]);
pdfplot(fcav);
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 13');
    print -dpdf ./bilder/C3_13.pdf 
end    
disp('Block = 16'), pause

%% Section 3.3.3 Rayleigh approximation for wave crest height
% Example 5: Transformed Rayleigh approximation of crest height  
% from spectral density
clf
xx = load('sea.dat');
x = xx;
x(:,2) = detrend(x(:,2));
SS = dat2spec2(x);
[sk, ku, me, si ] = spec2skew(SS);
gh = hermitetr([],[si sk ku me]);
Hs = 4*si;
r = (0:0.05:1.1*Hs)';
fac_h = trraylpdf(r,'Ac',gh);
fat_h = trraylpdf(r,'At',gh);
h = (0:0.05:1.7*Hs)';
facat_h = trraylpdf(h,'AcAt',gh);
pdfplot(fac_h);
hold on
pdfplot(fat_h,'--');
hold off
wafostamp([],'(ER)')
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 14');
    print -dpdf ./bilder/C3_14.pdf 
end    
disp('Block = 17'), pause

%%
%clf
TC = dat2tc(xx, me);
tc = tp2mm(TC);
Ac = tc(:,2);
At = -tc(:,1);
AcAt = Ac+At;
disp('Block = 18'), pause

%%
clf
Fac_h = [fac_h.x{1} cumtrapz(fac_h.x{1},fac_h.f)];
subplot(3,1,1)
Fac = plotedf(Ac,Fac_h);
hold on
plot(r,1-exp(-8*r.^2/Hs^2),'.')
axis([1. 2. 0.9 1])
title('Ac CDF')

Fat_h = [fat_h.x{1} cumtrapz(fat_h.x{1},fat_h.f)];
subplot(3,1,2)
Fat = plotedf(At,Fat_h);
hold on
plot(r,1-exp(-8*r.^2/Hs^2),'.')
axis([1. 2. 0.9 1])
title('At CDF')

Facat_h = [facat_h.x{1} cumtrapz(facat_h.x{1},facat_h.f)];
subplot(3,1,3)
Facat = plotedf(AcAt,Facat_h);
hold on
r2 = (0:0.05:2.1*Hs)';
plot(r2,1-exp(-2*r2.^2/Hs^2),'.')
axis([1.5 3.5 0.9 1])
title('At+Ac CDF')

wafostamp([],'(ER)'); hold off
if printing,
    annotation('TextBox',[0.05 0.94 0.14 0.05],...
        'FitBoxToText','on',...
        'String','Fig C3 15');
    print -dpdf ./bilder/C3_15.pdf 
end    
disp('Block = 19'), pause
disp('Elapsed time')
etime(clock,start)
