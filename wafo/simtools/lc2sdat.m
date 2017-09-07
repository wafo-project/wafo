function process=lc2sdat(lc,N,alpha)
%LC2SDAT Simulates process with given irregularity factor and crossing spectrum 
%
%  CALL: L = lc2sdat(lc,N,alpha);
%
%        L     = two column matrix with times and values of the simulated
%                process.
%        lc    = two column matrix with levels and crossing intensities,
%                i.e., level-crossing spectrum.. 
%        N     = number of sample points.
%        alpha = irregularity factor, 0<alpha<1, small  alpha  gives 
%                irregular process.
%
% Example
%  n = 10000;  
%  S = jonswap(7);
%  alpha = spec2char(S,'alpha');
%  xs  = spec2sdat(S,n);  
%  lc  = dat2lc(xs);
%  xs2 = lc2sdat(lc,n,alpha);
%  Se  = dat2spec(xs2);  
%  plotspec(S),hold on
%  plotspec(Se,'r'), hold off  
%  alpha2 = spec2char(Se,'alpha');
%  alpha - alpha2
%  lc2  = dat2lc(xs2);
%  figure(gcf+1)
%  subplot(211)
%  lcplot(lc2) 
%  subplot(212)
%  lcplot(lc)
  
%History
% revised pab Feb2004
% -added example  
% -changed order of input to conform with the other XXX2sdat routines  
% Revised by GL, July 10, 2000:
% Changed call, from  cross2tr  to  lc2tr
% Changed "Check the result" by removing reference to getrfc

% TODO % add a good example
  
%error(nargchk(1,3,nargin))
narginchk(1,3)
f = linspace(0,0.49999,1000);
rho_st = 2*sin(f*pi).^2-1;
tmp    = alpha*asin(sqrt((1+rho_st)/2));
tmp    = sin(tmp).^2;
a2     = (tmp-rho_st)./(1-tmp);
y      = min([a2+rho_st;1-a2]);
[maximum,maxidx]=max(y);

rho_st = rho_st(maxidx);
a2     = a2(maxidx);
a1     = 2*rho_st+a2-1;
r0     = 1;
r1     = -a1/(1+a2);
r2     = (a1^2-a2-a2^2)/(1+a2);
sigma2 = r0+a1*r1+a2*r2;

e      = randn(N,1)*sqrt(sigma2);
e(1:2) = [0;0];
L0     = randn(1,1);
L0     = [L0;r1*L0+sqrt(1-r2^2)*randn(1,1)];
%Simulate the process, starting in L0
L      = filter(1,[1 a1 a2],e,filter([1 a1 a2],1,L0));

epsilon = 1.01;
min_L   = min(L);
max_L   = max(L);
maxi    = max(abs([min_L max_L]))*epsilon;
mini    = -maxi;

u = linspace(mini,maxi,101)';
G = (1+erf(u/sqrt(2)))/2;
G = G.*(1-G);

x = linspace(0,r1,100)';
factor1  = 1./sqrt(1-x.^2);
factor2  = 1./(1+x);
integral = zeros(size(u));
for i=1:length(integral)
  y = factor1.*exp(-u(i)*u(i).*factor2);
  integral(i) = trapz(x,y);
end
G = G-integral/(2*pi);
G = G/max(G);

Z = ((u>=0)*2-1).*sqrt(-2*log(G));

sumcr   = trapz(lc(:,1),lc(:,2));
lc(:,2) = lc(:,2)/sumcr;
mcr     = trapz(lc(:,1),lc(:,1).*lc(:,2));
scr     = trapz(lc(:,1),lc(:,1).^2.*lc(:,2));
scr     = sqrt(scr-mcr^2);
g       = lc2tr(lc,mcr,scr);

f = [u u];
f(:,2) = tranproc(Z,fliplr(g));

process = tranproc(L,f);
process = [(1:length(process))' process];

%Check the result:
%save load.dat process -ascii
%[RFC TP]=getrfc;
%load rfc.dat
%param=[min(process(:,2)) max(process(:,2)) 100];
%fr=cc2fr(param,RFC);
%NT=fr2nt(fr);
%mu=pickdiag(NT);

%Check the result without reference to getrfc:
% LCe = dat2lc(process);
% max(lc(:,2))
% max(LCe(:,2))
% 
% clf
% plot(lc(:,1),lc(:,2)/max(lc(:,2)))
% hold on
% plot(LCe(:,1),LCe(:,2)/max(LCe(:,2)),'-.')
% title('Relative crossing intensity')

%% Plot made by the function funplot_4, JE 970707
%param = [min(process(:,2)) max(process(:,2)) 100];
%plot(lc(:,1),lc(:,2)/max(lc(:,2)))
%hold on
%plot(levels(param),mu/max(mu),'--')
%hold off
%title('Crossing intensity')  
%watstamp;

% Temporarily
%funplot_4(lc,param,mu);
