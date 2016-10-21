function [F_RFC] = refdemo2(demoNr,P,param,x0,s,lam)

% RFCDEMO2 Rainflow matrix for Switching Markov Chains of Turning Points.
%
% [F_RFC] = rfcdemo2;
% [F_RFC] = rfcdemo2(P);
% [F_RFC] = rfcdemo2(demoNr,P,param,x0,s,lam);
%
% Input:
% demoNr = 1:  Example 4.1  (Default: 1)
%          2:  Example 4.2
%          11: Another example
% P      = Transition matrix for regime process.  [rxr]
% param  = Parameter vector, [a b n], defines discretization.
% x0     = Center of ellipse [min Max]            [rx2]
% s      = Standard deviation (0<s<infinity)      [rx1]
% lam    = Orientation of ellipse (0<lam<infinity)[rx1]
%
% Output:
% F_RFC  = Calculated rainflow intensity          [nxn]
%
%   The demo presents calculation of the rainflow matrix for 
%   a Switching Markov Chains of Turning Points (SMCTP).
%   See Examples 4.1 and 4.2 in PhD-thesis:
%     P. Johannesson (1999): Rainflow Analysis of
%     Switching Markov Loads.
%     [http://www.maths.lth.se/matstat/staff/pj/PhD/]
%
% Example:
%   rfcdemo2;
%   rfcdemo2(1);
%   P=[0.98 0.02; 0.05 0.95]; 
%   rfcdemo2(1,P);
%   P=[0.98 0.02; 0.05 0.95]; x0=[-0.4 -0.3; 0.3 0.4]; param=[-1 1 32];
%   s=[0.1 0.15]'; lam=[0.7 1.5]'; 
%   rfcdemo2(1,P,param,x0,s,lam);
%
% See also  rfcdemo1, fatigue

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 18-May-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
% Updated by PJ 07-Jul-2005
%   warning off

warning_state=warning;  % Get current warning-settings
warning off

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(0,6,ni));


fprintf(1,'Welcome to demo 2 for "Rainflow Cycles for Switching Processes"!\n');
fprintf(1,'It demonstrates calculation of the rainflow matrix\n');
fprintf(1,'for a Switching Markov Chain of Turning Points.\n');
fprintf(1,'Type "help rfcdemo2" for further information\n\n');

Zstr='123456789';

% Delete all figure windows

slut=0;
h=gcf;
while ~slut
  h_prev=h;
  delete(h);
  h=gcf;
  if h==h_prev, slut=1; end
end
  
% Check input-arguments

if ni<1,        demoNr=1; end
if isempty(demoNr), demoNr=1; end
if ni<2,        P=[];     end
if ni<3,        param=[]; end
if ni<4,        x0=[];    end
if ni<5,        s=[];     end
if ni<6,        lam=[];   end

if demoNr == 1
  p1=0.10; p2=0.05;
  P0=[1-p1 p1; p2 1-p2];
  param0 = [-1 1 32];
  x00 = [-0.4 -0.3; 0.3 0.4];

  s0 = [0.15 0.15]';
  lam0 = [1.0 1.0]';
elseif demoNr == 2
  p1=0.10; p2=0.05;
  P0=[1-p1 p1; p2 1-p2];
  param0 = [-1 1 32];
  x00 = [-0.1 0.1; 0.0 0.0];

  s0 = [0.28 0.12]';
  lam0 = [0.5 2]';
elseif demoNr == 11
  p1=0.04; p2=0.01;
  P0=[1-p1 p1; p2 1-p2];
  param0 = [-1 1 32];
  x00 = [0.2 0.3; -0.1 -0.1];
  s0 = [0.18 0.1]';
  lam0 = [.5 2]';
end

if isempty(P),     P=P0;         end
if isempty(param), param=param0; end
if isempty(x0),    x0=x00;       end
if isempty(s),     s=s0;         end
if isempty(lam),   lam=lam0;     end

% Define min-Max and Max-min matrices

% Initiation

u=levels(param);    % Discretization levels
r=length(P);        % Number of regime states
n=param(3);         % Number of discretization levels

paramD = [1 n n];

F=cell(r,2);
for i=1:r
  [F{i,1},F{i,2}] = mktestmat(param,x0(i,:),s(i),lam(i));
end


% Calculate statistics: mean , std, proportion, mean time

%mX = m./sum(A')';
%sX = sqrt(s2./(1-A(:,2).^2));
statP = mc2stat(P);
tau = 1./(1-diag(P))';


% Parameters

fprintf('Parameters\n');
%fprintf(' z      a(z)      m(z)      s(z)\n');
%fprintf('%2d  %8g  %8g  %8g\n',[(1:r)' A(:,2) m sqrt(s2)]')
P

% Statistics
fprintf('Statistics\n');
fprintf('\n z    ro(z)    tau(z)\n');
fprintf('%2d  %8g  %8g\n',[(1:r)' statP' tau']')
fprintf('\n')

%fprintf('\n z    m_X(z)    s_X(z)     ro(z)    tau(z)\n');
%fprintf('%2d  %8g  %8g  %8g  %8g\n',[(1:r)' mX sX statP' tau']')

% Simulate a load

fprintf(1,' -- Simulating.\n');

T = 5000;   		% Number of TP

[xD,z] = smctpsim(P,F,T);
x = u(xD)';

I = 1:501;
plothmm(x(I),z(I),I-1,[1 r],'','Switching Process, X(k)',1);
ylabel('regime, Z(k)'),xlabel('time, k')

set(gcf,'Name','Switching MCTP')

drawnow, disp('Hit any key to continue.'); pause;

% Calculate the rainflow matrix

fprintf(1,' -- Calculating rainflow matrix.\n');

F_RFCsum = zeros(n,n);
for i=1:r
  [F_RFC1{i},mu_RFC1{i}] = smctp2rfm(1,F(i,1:2));
  F_RFCsum = F_RFCsum + statP(i)*F_RFC1{i};
end

figure
set(gcf,'Name','Calculated Rainflow matrices')
v = [param(1:2) param(1:2) 0 max(max(F_RFCsum))];
for i=1:r
  subplot(1,r,i),
  cmatplot(u,u,statP(i)*F_RFC1{i},1);
  grid,axis('square'),axis(v)
  title(['Regime ' Zstr(i)])
end

[F_RFC,mu_RFC] = smctp2rfm(P,F);

figure
set(gcf,'Name','Calculated Rainflow matrix')

v = [param(1:2) param(1:2) 0 max(max(F_RFCsum))];
subplot(1,2,1)
cmatplot(u,u,F_RFC,1),grid
axis('square'),axis(v)
title('Mixed RFM')
subplot(1,2,2)
cmatplot(u,u,F_RFCsum,1),grid
axis('square'),axis(v)
title('Sum of RFMs')

drawnow, disp('Hit any key to continue.'); pause;

figure

% Calculate Cycles and crossings

fprintf(1,' -- Cycles and crossings for simulated load.\n');

TP = dat2tp([(1:T)' x]);
RFC = tp2rfc(TP);
plotcc(RFC)
cross = cc2lc(RFC,1);   % Only upcrossings

% Observed Rainflow Matrix

TPD = dat2tp([(1:T)' xD]);
RFCD = tp2rfc(TPD);
fr = cc2cmat(paramD,RFCD);
F_RFCobs = fr;
%[F_RFCobs_smooth,h] = smthcmat(F_RFCobs);
[F_RFCobs_smooth,h] = smoothcmat(F_RFCobs,1,0.7);

set(gcf,'Name','Observed Rainflow Matrix')
subplot(1,2,1), cmatplot(u,u,F_RFCobs/T*2,1),grid
axis('square'),axis(v)
title('Observed RFM')
subplot(1,2,2), cmatplot(u,u,F_RFCobs_smooth/T*2,1), grid
axis('square'),axis(v)
title('Smoothed observed RFM')

drawnow, disp('Hit any key to continue.'); pause;


% Check crossing intensity

% Check Damage matrix

figure
set(gcf,'Name','Damage')
beta0 = 3;
Dmat = cmat2dmat(param,F_RFC,beta0);
Dmatsum = cmat2dmat(param,F_RFCsum,beta0);
subplot(1,2,1),cmatplot(u,u,Dmat),grid
axis('square')

v=axis;
title('Mixed RFM')
subplot(1,2,2),cmatplot(u,u,Dmatsum),grid
axis('square')

axis(v);
title('Sum of RFMs')


drawnow, disp('Hit any key to continue.'); pause;


% Compare crossing spectrum with simulated 

figure
set(gcf,'Name','Crossing Intensity')

subplot(1,2,1)
plot(u,diag(mu_RFC),'g--'), hold on
%[XX,YY]=stairs(u,diag(mu_RFC));
%plot(XX,YY,'g'),hold on
%stairs(cross(:,1),cross(:,2)/T*2), hold off
plot(cross(:,1),cross(:,2)/T*2), hold off
xlabel('level, u'); ylabel('mu(u)');
v = axis; axis([param(1:2) 0 v(4)]);

subplot(1,2,2)
plot(u,diag(mu_RFC),'g'), hold on
plot(u,diag(mu_RFC1{1})*statP(1),'r-.')
plot(u,diag(mu_RFC1{2})*statP(2),'r-.')
plot(u,diag(mu_RFC1{1})*statP(1)+diag(mu_RFC1{2})*statP(2),'y--'), hold off
xlabel('level, u'); ylabel('mu(u)');
v = axis; axis([param(1:2) 0 v(4)]);

fprintf(1,'\nThank you for running this demo!\n');

warning(warning_state);


