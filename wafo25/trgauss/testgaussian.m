function test2 = testgaussian(S,Np,test0,def,opt)
%TESTGAUSSIAN Test if a stochastic process is Gaussian.
%
% CALL:  test1 = testgaussian(S,[Np,Ns],test0,def,options);
%
% test1,
%    test0 = simulated and observed value of e(g)=int (g(u)-u)^2 du,
%            respectively, where int limits is given by OPTIONS.PARAM. 
%         
%      S   = spectral density structure 
%
%       Np = # of points simulated
%       Ns = # of independent simulations (default  100)
%
%   def    = 'nonlinear' : transform based on smoothed crossing intensity (default)
%            'mnonlinear': transform based on smoothed marginal distribution
%  options = options structure defining how the estimation of the
%            transformation is done. (default troptset('dat2tr'))
%
% TESTGAUSSIAN simulates  e(g(u)-u) = int (g(u)-u)^2 du  for Gaussian processes 
% given the spectral density, S. The result is plotted if test0 is given.
% This is useful for testing if the process X(t) is Gaussian.
% If 95% of TEST1 is less than TEST0 then X(t) is not Gaussian at a 5% level.
%
% Example:
%    Hm0 = 7;
%    S0 = jonswap([],Hm0); g=ochitr([],[Hm0/4]); S=S0;
%    S.tr=g;S.tr(:,2)=g(:,2)*Hm0/4;
%    xs = spec2sdat(S,2^13);
%    [g0 t0] = dat2tr(xs);
%    t1 = testgaussian(S0,[2^13 50],t0); 
%
% See also  cov2sdat, dat2tr, troptset

%Tested on: Matlab 5.3, 5.2, 5.1
% History:
% revised pab
% -changed name from mctrtest to testgaussianity
% revised pab jan2005
% changed order of input so that chapter1.m works
% revised pab Feb2004
% -changed h1line  
% revised pab 29.12.2000
% - added options and def to input due to changed calling syntax for dat2tr.
% - updated the help header
% revised by pab 12.11.99
%    fixed a bug, string input to dat2tr 'nonlinear' 
% revised by pab 12.10.99
%    updated help header 
% revised by pab 11.08.99
% changed name from mctest to mctrtest
% by pab 11.11.98

maxsize = 200000; % must divide the computations due to limited memory


error(nargchk(2,5,nargin))
if nargin<4||isempty(def),
  def = 'nonlinear';
end
if nargin<5||isempty(opt),
  opt = troptset('dat2tr');
end
opt = troptset(opt,'multip',1);

if nargin<3||isempty(test0)
  plotflag=0;
else 
  plotflag=1;
end

Ns=100;
nnp=length(Np);
if nnp>=2 
  Ns=Np(2);
end
Np=Np(1);
if Ns>50
  disp('  ... be patient this may take a while')
end

test1 = [];
rep = floor(Np*Ns/maxsize)+1;

Nstep = floor(Ns/rep);

R = spec2cov(S);

for ix=1:rep,
  xs = cov2sdat(R,[Np Nstep]);
  [g, tmp] = dat2tr(xs,def,opt);
  test1 = [test1; tmp(:)];
  disp(['finished ' num2str(ix) ' of ' num2str(rep)] )
end
if rep>1,
  xs = cov2sdat(R,[Np rem(Ns,rep)]);
  [g tmp] = dat2tr(xs,def,opt);
  test1 = [test1; tmp(:)];
end

if (nargout>0 || plotflag==0),
  test2=test1;
end


if plotflag 
  plot(test1,'o'),hold on
  if 1 
    plot([1 Ns],test0*[1 1],'--'),
  end
  hold off
  ylabel('e(g(u)-u)')
  xlabel('Simulation number')
end
