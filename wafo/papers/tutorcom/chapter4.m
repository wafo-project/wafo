%% CHAPTER4  Demonstrates exact distributions of wave characteristics
%
% Chapter4 contains the commands used in Chapter4 in the tutorial.
% 
% Some of the commands are edited for fast computation. 
% Each set of commands is followed by a 'pause' command.
% Set pstate='on' to activate the pause option

% Tested on Matlab 5.3, 7.10
% History
% Revised by Georg Lindgren march 2011 for Tutorial 2.5 and 
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
%pstate = 'on';
speed = 'fast';
%speed = 'slow'
pause(pstate)

%% Section 4.2 Marginal distributions of wave characteistics
%% Section 4.2.1 Crest period, crest length, and crest height

% Example 6 Torsethaugen waves
clf
S1 = torsethaugen([],[6 8],1);
D1 = spreading(101,'cos',pi/2,[15],[],0);
D12 = spreading(101,'cos',0,[15],S1.w,1);
SD1 = mkdspec(S1,D1);
SD12 = mkdspec(S1,D12);
disp('Block = 1'), pause

%% 6a. Crest period
clf
tic 
f_tc_4 = spec2tpdf(S1,[],'Tc',[0 12 61],[],4);
f_tc_1 = spec2tpdf(S1,[],'Tc',[0 12 61],[],-1);
toc
pdfplot(f_tc_4,'-.'), hold on
pdfplot(f_tc_1), hold off
simpson(f_tc_4.x{1},f_tc_4.f)
simpson(f_tc_1.x{1},f_tc_1.f)
wafostamp([],'(ER)')
disp('Block = 2'), pause

%% Setting of rind options

if strncmpi(speed,'slow',1)
  opt1 = rindoptset('speed',5,'method',3);
  opt2 = rindoptset('speed',5,'nit',5,'method',0);
else
  % fast
  opt1 = rindoptset('speed',7,'method',3);
  opt2 = rindoptset('speed',7,'nit',2,'method',0);
end

clf
if strncmpi(speed,'slow',1)
  NITa = 5;
else
  disp('NIT=5 may take time, running with NIT=3 in the following')
  NITa = 3;
end

%% 6b. Crest length

% NIT = -1 gives fast computations - you should compare with slower 
% alternatives, NIT = 3 or 5 that give upper bounds
%f_Lc = spec2tpdf2(S1,[],'Lc',[0 200 81],opt1); % Faster and more accurate

clf
f_Lc = spec2tpdf(S1,[],'Lc',[0 125 251],[],-1);
pdfplot(f_Lc,'-.'), hold on
wafostamp([],'(ER)')
disp('Block = 3'), pause

f_Lc_1 = spec2tpdf(S1,[],'Lc',[0 125 251],1.5,-1);
%f_Lc_1 = spec2tpdf2(S1,[],'Lc',[0 200 81],1.5,opt1);
hold on
pdfplot(f_Lc_1), hold off
wafostamp([],'(ER)')

disp('Block = 4'), pause
%% 
simpson(f_Lc.x{1},f_Lc.f)
simpson(f_Lc_1.x{1},f_Lc_1.f)
      
disp('Block = 5'), pause

%% 6c. Crest height
% Crest height cdf from spec2acdf compared to simulated data and the 
% Rayleigh approximation with  Hs=6

clf
Hs=6;
r = (0:0.06:Hs)'; %Compute cdf-values at r
f_Ac_s1_T = spec2acdf(S1,[],'Tc',[0 12 61],r,-1); 
hold on
T = spec2sdat(S1,[40000,100],0.01);
[Steep,Height,AcT] = dat2steep(T);
plotedf(AcT,'-.')
if strncmpi(speed,'slow',1), 
    Lev=[0 120 241];
else
    Lev=[0 120 61];
end
F_Ac_s1_L = spec2acdf(S1,[],'Lc',Lev,r,-1);
L = spec2sdat(spec2spec(S1,'k1d'),[40000 100],0.1);
[Steep,Height,AcL] = dat2steep(L);
plotedf(AcL,'-.'), hold off

disp('Block = 6'), pause

%% 6d. Directional spreading

% pdf of Lc

figure(1)
clf
if strncmpi(speed,'slow',1), 
    Lev2=[0 200 401];
else
    Lev2=[0 200 201];
end
tic
f_Lc_d1 = spec2tpdf(rotspec(SD1,pi/2),[],'Lc',Lev2,[],-1);
f_Lc_d12 = spec2tpdf(SD12,[],'Lc',Lev2,[],-1);
% f_Lc_d1 = spec2tpdf2(rotspec(SD1,pi/2),[],'Lc',[0 200 81],opt1);
% f_Lc_d12 = spec2tpdf2(SD12,[],'Lc',[0 200 81],opt1);
toc
pdfplot(f_Lc_d1,'-.'), hold on
pdfplot(f_Lc_d12),     hold off
wafostamp([],'(ER)')

disp('Block = 7'), pause

figure(2)
dx = f_Lc.x{1}(2)-f_Lc.x{1}(1);  
dx1 = f_Lc_d1.x{1}(2)-f_Lc_d1.x{1}(1);  
dx12 = f_Lc_d12.x{1}(2)-f_Lc_d12.x{1}(1);  
plot(f_Lc.x{1},cumsum(f_Lc.f)*dx), hold on
plot(f_Lc_d1.x{1},cumsum(f_Lc_d1.f)*dx1,'-.')
plot(f_Lc_d12.x{1},cumsum(f_Lc_d12.f)*dx12,'--')
hold off
wafostamp([],'(ER)')
disp('Block = 8'), pause

%% Section 4.2.2 Numerical accuracy
%%

figure(1)
clf
opt1 = rindoptset('speed',5,'method',3);
SD1r = rotspec(SD1,pi/2);
if strncmpi(speed,'slow',1)
  f_Lc_d1_5 = spec2tpdf(SD1r,[], 'Lc',[0 200 201],[],5);
  pdfplot(f_Lc_d1_5),     hold on
else
  % fast
  disp('Run the following example only if you want a check on computing time')
  disp('Edit the command file and remove %')
end
f_Lc_d1_3 = spec2tpdf(SD1r,[],'Lc',[0 200 201],[],3);
f_Lc_d1_2 = spec2tpdf(SD1r,[],'Lc',[0 200 201],[],2);
f_Lc_d1_0 = spec2tpdf(SD1r,[],'Lc',[0 200 201],[],0);
f_Lc_d1_neg = spec2tpdf(SD1r,[],'Lc',[0 200 201],[],-1);
%f_Lc_d1_n4 = spec2tpdf2(SD1r,[],'Lc',[0 400 161],opt1);

pdfplot(f_Lc_d1_3), hold on
pdfplot(f_Lc_d1_2)
pdfplot(f_Lc_d1_0)
pdfplot(f_Lc_d1_neg)
%pdfplot(f_Lc_d1_n4)
hold off
%simpson(f_Lc_d1_n4.x{1},f_Lc_d1_n4.f)

disp('Block = 9'), pause

%% Section 4.2.3 Wave period and wave length
% The wave period (length) is the sum of crest period and the following 
% trough period (length) and is complicated to compute
%% Example 7: Crest period and high crest waves
clf
tic
xx = load('sea.dat');
x = xx;
x(:,2) = detrend(x(:,2));
SS = dat2spec(x);
si = sqrt(spec2mom(SS,1));
SS.tr = dat2tr(x);
Hs = 4*si

%% 7a. Crest period (again)
method = 0;
rate = 2;
[S, H, Ac, At, Tcf, Tcb, z_ind, yn] = dat2steep(x,rate,method);
Tc = Tcf+Tcb;
t = linspace(0.01,8,200);
f_tc1emp = kde(Tc,{'L2',0},t);
pdfplot(f_tc1emp)
hold on
if strncmpi(speed,'slow',1)
  NIT = 5;
else
  NIT = 2;
end
f_tc1 = spec2tpdf(SS,[],'Tc',[0 8 81],0,NIT);
simpson(f_tc1.x{1},f_tc1.f)
pdfplot(f_tc1,'-.')
hold off
wafostamp([],'(ER)')
toc
disp('Block = 10'), pause

%% Crest period for high crests
clf
if strncmpi(speed,'slow',1)
  NIT = 5;
else
  NIT = 2;
end
tic
f_tc2 = spec2tpdf(SS,[],'Tc',[0 8 81],Hs/2,NIT);
toc
Pemp = sum(Ac>Hs/2)/sum(Ac>0)
simpson(f_tc2.x{1},f_tc2.f)
index = find(Ac>Hs/2);
f_tc2emp = kde(Tc(index),{'L2',0},t);
f_tc2emp.f = Pemp*f_tc2emp.f;
pdfplot(f_tc2emp)
hold on
pdfplot(f_tc2,'-.')
hold off
wafostamp([],'(ER)')
toc
disp('Block = 11'), pause

%% 7b. Wave period for high crest waves 
% This is moderately hard
      clf
      tic 
      f_tu1_n = spec2tccpdf(SS,[],'t>',[0 12 61],[Hs/2],[0],-1);
      toc
      simpson(f_tu1_n.x{1},f_tu1_n.f)
      f_tu1_3 = spec2tccpdf(SS,[],'t>',[0 12 61],[Hs/2],[0],3,5);
%      f_tu1_5 = spec2tccpdf(SS,[],'t>',[0 12 61],[Hs/2],[0],5,5);
      simpson(f_tu1_3.x{1},f_tu1_3.f)
      pdfplot(f_tu1_n,'-.')
      hold on
      pdfplot(f_tu1_3)
      hold off
      toc
disp('Block = 12'), pause

%% 7c. Wave period for high-crest, deep-trough waves
% This is a test for accuracy
      clf
      [TC tc_ind v_ind] = dat2tc(yn,[],'dw');
      N = length(tc_ind);
      t_ind = tc_ind(1:2:N);
      c_ind = tc_ind(2:2:N);
      Pemp = sum(yn(t_ind,2)<-Hs/2 & ...
          yn(c_ind,2)>Hs/2)/length(t_ind)
      ind = find(yn(t_ind,2)<-Hs/2 & yn(c_ind,2)>Hs/2);
      spwaveplot(yn,ind(2:4))
      wafostamp([],'(ER)')
disp('Block = 13'), pause

%% Upcrossing period Tu for high crest, deep trough waves 

clf
Tu = yn(v_ind(1+2*ind),1)-yn(v_ind(1+2*(ind-1)),1);
t = linspace(0.01,14,200);
f_tu2_emp = kde(Tu,{'kernel' 'epan','L2',0},t);
f_tu2_emp.f = Pemp*f_tu2_emp.f;
pdfplot(f_tu2_emp,'-.')
wafostamp([],'(ER)')
disp('Block = 14'), pause

tic 
f_tu2_n = spec2tccpdf(SS,[],'t>',[0 12 61],[Hs/2],[Hs/2],-1);
toc
simpson(f_tu2_n.x{1},f_tu2_n.f)
hold on
pdfplot(f_tu2_n)
hold off
wafostamp([],'(ER)')
disp('Block = 15'), pause

%%
disp('The rest of this chapter deals with joint densities.')
disp('Some calculations may take some time.') 
disp('You could experiment with other NIT.')
%return

%% Section 4.3 Joint density of crest period and crest height
%% Section 4.3.1. Some preliminary analysis of the data
% Example 8. Gullfaks
clf
tic
yy = load('gfaksr89.dat');
SS = dat2spec(yy);
si = sqrt(spec2mom(SS,1));
SS.tr = dat2tr(yy);
Hs = 4*si
v = gaus2dat([0 0],SS.tr);
v = v(2)
toc
disp('Block = 16'), pause

%%
clf
tic
[TC, tc_ind, v_ind] = dat2tc(yy,v,'dw');
N       = length(tc_ind);
t_ind   = tc_ind(1:2:N);
c_ind   = tc_ind(2:2:N);
v_ind_d = v_ind(1:2:N+1);
v_ind_u = v_ind(2:2:N+1);
T_d     = ecross(yy(:,1),yy(:,2),v_ind_d,v);
T_u     = ecross(yy(:,1),yy(:,2),v_ind_u,v);

Tc = T_d(2:end)-T_u(1:end);
Tt = T_u(1:end)-T_d(1:end-1);
Tcf = yy(c_ind,1)-T_u;
Ac = yy(c_ind,2)-v;
At = v-yy(t_ind,2);
toc
disp('Block = 17'), pause

%%
clf
tic
t = linspace(0.01,15,200);
kopt3 = kdeoptset('hs',0.25,'L2',0); 
ftc1 = kde(Tc,kopt3,t);
ftt1 = kde(Tt,kopt3,t);
pdfplot(ftt1,'k')
hold on
pdfplot(ftc1,'k-.')
f_tc4 = spec2tpdf(SS,[],'Tc',[0 12 81],0,4,5);
%f_tc2 = spec2tpdf(SS,[],'Tc',[0 12 81],0,2,5);
f_tcn = spec2tpdf(SS,[],'Tc',[0 12 81],0,-1);
pdfplot(f_tcn,'b')
hold off
legend('kde(Tt)','kde(Tc)','f_{tc}')
wafostamp([],'(ER)')
toc
disp('Block = 18'), pause

%% Section 4.3.2 Joint distribution of crest period and height
%% Example 9. Joint characteristics of a half wave:
%% position and height of a crest for a wave with given period
clf
tic
ind = find(4.4<Tc & Tc<4.6);
f_AcTcf = kde([Tcf(ind) Ac(ind)],{'L2',[1 .5]});
pdfplot(f_AcTcf)
hold on
plot(Tcf(ind), Ac(ind),'.');
wafostamp([],'(ER)')
toc
disp('Block = 19'), pause

%%
clf
tic
opt1 = rindoptset('speed',5,'method',3);
opt2 = rindoptset('speed',5,'nit',2,'method',0);

f_tcfac1 = spec2thpdf(SS,[],'TcfAc',[4.5 4.5 46],[0:0.25:8],opt1);
f_tcfac2 = spec2thpdf(SS,[],'TcfAc',[4.5 4.5 46],[0:0.25:8],opt2);

pdfplot(f_tcfac1,'-.')
hold on
pdfplot(f_tcfac2)
plot(Tcf(ind), Ac(ind),'.'); hold off

simpson(f_tcfac1.x{1},simpson(f_tcfac1.x{2},f_tcfac1.f,1))
simpson(f_tcfac2.x{1},simpson(f_tcfac2.x{2},f_tcfac2.f,1))
f_tcf6=spec2tpdf(SS,[],'Tc',[4.5 4.5 46],[0:0.25:8],6);
f_tcf6.f(46)
toc
wafostamp([],'(ER)')
disp('Block = 20'), pause

%% 
if strncmpi(speed,'slow',1)
clf
tic
f_tcac_s = spec2thpdf(SS,[],'TcAc',[0 12 81],[Hs/2:0.1:2*Hs],opt1);
toc
disp('Block = 21'), pause

clf
tic
mom = spec2mom(SS,4,[],0);
t = f_tcac_s.x{1};
h = f_tcac_s.x{2};
flh_g = lh83pdf(t',h',[mom(1),mom(2),mom(3)],SS.tr);
clf
ind=find(Ac>Hs/2);
plot(Tc(ind), Ac(ind),'.');
hold on
pdfplot(flh_g,'k-.')
pdfplot(f_tcac_s)
toc
wafostamp([],'(ER)')
disp('Block = 22'), pause
end
%%
clf
%      f_tcac = spec2thpdf(SS,[],'TcAc',[0 12 81],[0:0.2:8],opt1);
%      pdfplot(f_tcac)
disp('Block = 23'), pause

%% Section 4.3.3 Joint density of crest and trough height
%% Section 4.3.4 Min-to-max distributions – Markov method
%% Example 10. (min-max problems with Gullfaks data)
%% Joint density of maximum and the following minimum
clf
tic
opt2 = rindoptset('speed',5,'nit',2,'method',0);
tp = dat2tp(yy);
Mm = fliplr(tp2mm(tp));
fmm = kde(Mm);
f_mM = spec2mmtpdf(SS,[],'mm',[],[-7 7 51],opt2);

pdfplot(f_mM,'-.')
hold on
pdfplot(fmm,'k-')
hold off
wafostamp([],'(ER)')
toc
disp('Block = 24'), pause

%% The joint density of ”still water separated”  maxima and minima.
%% Example 11. crest-trough dirstribution from max-min transitions
clf
tic
opt2 = rindoptset('speed',5,'nit',2,'method',0);
ind = find(Mm(:,1)>v & Mm(:,2)<v);
Mmv = abs(Mm(ind,:)-v);
fmmv = kde(Mmv);
f_vmm = spec2mmtpdf(SS,[],'vmm',[],[-7 7 51],opt2);
clf
pdfplot(fmmv,'k-')
hold on
pdfplot(f_vmm,'-.')
hold off
wafostamp([],'(ER)')
toc
disp('Block = 25'), pause


%%
clf
tic
facat = kde([Ac At]);
f_acat = spec2mmtpdf(SS,[],'AcAt',[],[-7 7 51],opt2);
clf
pdfplot(f_acat,'-.')
hold on
pdfplot(facat,'k-')
hold off
wafostamp([],'(ER)')
toc
disp('Block = 26'), pause

disp('Elapsed time')
etime(clock,start)
