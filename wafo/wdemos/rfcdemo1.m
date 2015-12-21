function [F_RFC] = rfcdemo1(demoNr,P,A,m,s2,param)
% RFCDEMO1  Demo for switching AR(1)-processes. 
%
% [F_RFC] = rfcdemo1;
% [F_RFC] = rfcdemo1(demoNr);
% [F_RFC] = rfcdemo1(demoNr,P);
% [F_RFC] = rfcdemo1(demoNr,P,A,m,s2,param);
%
% Input:
% demoNr = 1 / 2 / 3 / 0
%          Default: 1. 0 = Define your own process.
% P      = Transition matrix for regime process.  [rxr]
% A      = AR coefficients.                       [rx2]
% m      = Mean koefficients.                     [rx1]
% s2     = Innovation variances.                  [rx1]
% param  = Parameter vector, [a b n], defines discretization.
%
% Output:
% F_RFC  = Calculated rainflow intensity          [nxn]
%
%   Demo for switching AR(1)-processes, by using 
%   switching Markov chains. 
%   See Examples 3.1-3.3 in PhD-thesis:
%     P. Johannesson (1999): Rainflow Analysis of
%     Switching Markov Loads.
%     [http://www.maths.lth.se/matstat/staff/pj/PhD/]
%
% Example:
%   rfcdemo1;
%   rfcdemo1(2);
%   P=[0.98 0.02; 0.05 0.95]; rfcdemo1(1,P);
%
% See also  rfcdemo2, fatigue

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Pär Johannesson) 1997
% Updated by PJ 19-May-2000
%   updated for WAFO
% Updated by PJ 07-Jul-2005
%   warning off

warning_state=warning;  % Get current warning-settings
warning off

fprintf(1,'Welcome to demo 1 for "Rainflow Cycles for Switching Processes"!\n');
fprintf(1,'It demonstrates calculation of the rainflow intensity\n');
fprintf(1,'for a switching AR(1)-process.\n');
fprintf(1,'Type "help rfcdemo1" for further information\n\n');

% Check input arguments

ni = nargin;
%no = nargout;
error(nargchk(0,6,ni));

Zstr='123456789';

% Add path to rfcdemo1
demoDir = [waforoot filesep 'wdemos' filesep 'rfcdemo1'];
addpath(demoDir);

% Delete all figure windows

slut=0;
h=gcf;
while ~slut
  h_prev=h;
  delete(h);
  h=gcf;
  if h==h_prev, slut=1; end
end
  

if ni == 0, demoNr=1; end

% Define switching AR(1)-process

if demoNr == 1
  p1=0.02; p2=0.01;
  P0 = [1-p1 p1; p2 1-p2];
  A = [1 -0.5; 1 0.5];
  m = 2*[-0.5 1.5]';
  s2 = [1 1]';
  param = [-8 8 64];
elseif demoNr == 2
  p1=0.01; p2=0.02;
  P0 = [1-p1 p1; p2 1-p2];
  A = [1 -0.5; 1 0.5];
  m = [0 0]';
  s2 = [4 1]';
  param = [-10 10 64];
elseif demoNr == 3
  p12=0.02; p13=0.005;
  p21=0.01; p23=0.01;
  p31=0.005; p32=0.02;
  P0 = [1-p12-p13 p12 p13; p21 1-p21-p23 p23; p31 p32 1-p31-p32];
  A = [1 -0.5; 1 -0.3; 1 0.5];
  m = 2*[-0.5 0 1.5]';
  s2 = [1 1 1.44]';
  param = [-8 8 64];
elseif demoNr == 0
  if ni < 6
    error('demoNr=0: Too few input arguments.');
  end
end

if ~exist('P','var')
  P=P0;
end

% Initiation

u=levels(param);    % Discretization levels
delta = u(2)-u(1);  % Discretization step
r=length(P);        % Number of regime states
n=param(3);         % Number of discretization levels
C=ones(r,1);        % 

% Calculate statistics: mean , std, proportion, mean time

mX = m./sum(A')';
sX = sqrt(s2./(1-A(:,2).^2));
statP = mc2stat(P);
tau = 1./(1-diag(P))';

% Parameters

fprintf('Parameters\n');
fprintf(' z      a(z)      m(z)      s(z)\n');
fprintf('%2d  %8g  %8g  %8g\n',[(1:r)' A(:,2) m sqrt(s2)]')
P
fprintf('\n')

% Statistics
fprintf('Statistics\n');
fprintf('\n z    m_X(z)    s_X(z)     ro(z)    tau(z)\n');
fprintf('%2d  %8g  %8g  %8g  %8g\n',[(1:r)' mX sX statP' tau']')
fprintf('\n')

% Calculate spectral densities

fprintf(1,' -- Calculate spectral densities.\n');

set(gcf,'Name','Power spectra')

for i=1:r
  Si = armaspec(1,A(i,:),s2(i,:),100);
  subplot(r,1,i)
  plot(Si(:,1),Si(:,2)), grid
  ylabel(['S(f), regime ' Zstr(i)])
end
xlabel('frequency, f')

drawnow, disp('Hit any key to continue.'); pause;

% Simulate MARX(2,1)-process

fprintf(1,' -- Simulating.\n');

T=10000;   % Length of simulation
Tinit=500; % Length of initiation (to remove transients)

[x,z,e] = sarmasim(C,A,m,s2,P,T,Tinit);

figure
set(gcf,'Name','Switching AR(1)-process')

I = 1:501;
hmmplot(x(I),z(I),I-1,[1 r],'','Switching Process, X(k)',1);
ylabel('regime, Z(k)'),xlabel('time, k')

drawnow, disp('Hit any key to continue.'); pause;

%
% Discretize AR(1)-processes
% Approximate with Markov chain
%

ud = u; ud(1)=u(1)-0.25*(u(n)-u(1)); ud(n)=u(n)+0.25*(u(n)-u(1));

figure
set(gcf,'Name','Markov chain for Discretized AR(1)-processes')

if demoNr == 1 || demoNr ==2 || demoNr ==3
  fprintf(1,' -- Loading discretized AR-processes.\n');

  eval(['load rfcdem1' Zstr(demoNr)]); % load uu,Q1,Q2,Q3
  for i=1:r
    subplot(1,r,i)
    eval(['Q{i}=Q' num2str(i) ';']);
    cmatplot(uu,uu,Q{i},3), axis square
    title(['Regime ' Zstr(i)])
  end;
  drawnow
else

  fprintf(1,' -- Discretizing AR-processes. (This can take some time!)\n');
  for i=1:r
    fprintf(1,' -- Discretizing process %d.\n',i);
    a1i = A(i,2);
    si = sqrt(s2(i));
    mi = mX(i);
    [Qi,uu] = discar(ud,a1i,mi,si);
    Q{i}=Qi;
    %eval(['Q{i}' num2str(i) '=Qi;']);
    subplot(r,1,i)
    cmatplot(uu,uu,Qi,3), axis square
    title(['Regime ' Zstr(i)])
    drawnow
  end
end

drawnow, disp('Hit any key to continue.'); pause;

%
% Compute RFC counting distribution for SMC
%

fprintf(1,' -- Calculating rainflow intensity.\n');

Qstr='';
for i=1:r
  Qstr = [Qstr ',Q' Zstr(i)];
end

[F_RFC,mu_RFC] = smc2rfm(P,Q,[2 delta]);    % definition 2
%eval( ['[F_RFC,mu_RFC] = smc2rfc(P,[2 delta]' Qstr ');'] );    % definition 2

for i=1:r
  %eval( ['[F_RFC' Zstr(i) ',mu_RFC' Zstr(i) '] = mc2rfc(Q' Zstr(i) ',[2 delta]);'] );     % Regime i
  [F_RFC0{i},mu_RFC0{i}] = mc2rfm(Q{i},[2 delta]);     % Regime i
end

figure

cmatplot(u,u,F_RFC,2);
set(gcf,'Name','Rainflow intensity for switching AR(1)-process')
title('3D-plot')

drawnow, disp('Hit any key to continue.'); pause;

figure

% Calculate Cycles and crossings

fprintf(1,' -- Cycles and crossings for simulated load.\n');

TP = dat2tp([(1:T)' x]);
RFC = tp2rfc(TP);
cross = cc2lc(RFC);

cocc(param,RFC,F_RFC)
set(gcf,'Name','Rainflow intensity for switching AR(1)-process')
title('Iso-lines compared with observed cycles')
xlabel('min'), ylabel('Max')

drawnow, disp('Hit any key to continue.'); pause;

% Plot Crossing intensity

figure
set(gcf,'Name','Crossing intensity')

subplot(1,2,1)
plot(u,diag(mu_RFC),'g--'), hold on
stairs(cross(:,1),cross(:,2)/T), hold off
title('Calculated vs. observed')
xlabel('level, u'); ylabel('Crossing intensity, mu(u)');

subplot(1,2,2)
plot(u,diag(mu_RFC),'g'), hold on
cross_sum=zeros(size(diag(mu_RFC)));
for i=1:r
  crossi=diag(mu_RFC0{i});
  plot(u,statP(i)*crossi,'r-.')
  cross_sum=cross_sum+statP(i)*crossi;
end
plot(u,cross_sum,'y--'), hold off
title('Calculated vs. individual')
xlabel('level, u'); ylabel('Crossing intensity, mu(u)');

fprintf(1,'\nThank you for running this demo!\n');

% Remove path to rfcdemo1
demoDir = [waforoot filesep 'wdemos' filesep 'rfcdemo1'];
rmpath(demoDir);

warning(warning_state);

