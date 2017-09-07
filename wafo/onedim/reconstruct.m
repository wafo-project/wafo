function [y,g,g2,test,tobs,mu1o, mu1oStd]=reconstruct(x,inds,Nsim,L,def,varargin)
%RECONSTRUCT reconstruct the spurious/missing points of timeseries
%
% CALL: [y,g,g2,test,tobs,mu1o,mu1oStd]=reconstruct(x,inds,Nsim,L,def,options);
%   
%       y   = reconstructed signal
%      g,g2 = smoothed and empirical transformation, respectively
% test,tobs = test observator int(g(u)-u)^2 du and int(g_new(u)-g_old(u))^2 du,
%             respectively, where int limits is given by param in lc2tr. 
%             Test is a measure of departure from the Gaussian model for 
%             the data. Tobs is a measure of the convergence of the 
%             estimation of g.
%      mu1o = expected surface elevation of the Gaussian model process.
%   mu1oStd = standarddeviation of mu1o.
%
%       x   = 2 column timeseries 
%             first column sampling times [sec]
%             second column surface elevation [m]
%      inds = indices to spurious points of x
%      Nsim = the maximum # of iterations before we stop
%
%         L = lag size of the Parzen window function. 
%             If no value is given the lag size is set to
%             be the lag where the auto correlation is less than 
%             2 standard deviations. (maximum 200) 
%       def = 'nonlinear' : transform based on smoothed crossing intensity (default)
%             'mnonlinear': transform based on smoothed marginal distribution
%             'linear'    : identity.
%   options = options structure defining how the estimation of g is
%             done, see troptset.
%   
% In order to reconstruct the data a transformed Gaussian random process is
% used for modelling and simulation of the missing/removed data conditioned
% on the other known observations.
%
% Estimates of standarddeviations of y is obtained by a call to tranproc
%        Std = tranproc(mu1o+/-mu1oStd,fliplr(g));
%
% See also  troptset, findoutliers, cov2csdat, dat2cov, dat2tr, detrendma

% Reference
% Brodtkorb, P, Myrhaug, D, and Rue, H (2001)
% "Joint distribution of wave height and wave crest velocity from
% reconstructed data with application to ringing"
% Int. Journal of Offshore and Polar Engineering, Vol 11, No. 1, pp 23--32 
%
% Brodtkorb, P, Myrhaug, D, and Rue, H (1999)
% "Joint distribution of wave height and wave crest velocity from
% reconstructed data"
% in Proceedings of 9th ISOPE Conference, Vol III, pp 66-73

% tested on: Matlab 5.3, 5.1
% History:
% revised pab 18.04.2001 
% - updated help header by adding def to function call
% - fixed a bug: param was missing
% revised pab 29.12.2000
% - replaced csm1,...., param, with a options structure
% - monitor replaced with plotflag==2
% revised pab 09.10.2000
% - updated call to dat2cov
% revised pab 22.05.2000
% - found a bug concerning cvmax and indr
% - updated call to waveplot
% revised pab 27.01.2000
% - added L,csm1,csm2,param, monitor to the input arguments 
% revised pab 17.12.1999
% -added reference
% revised pab 12.10.1999
%   updated arg. list to dat2tr
% last modified by Per A. Brodtkorb 01.10.98 

opt = troptset('dat2tr');
if nargin==1 && nargout <= 1 && isequal(x,'defaults')
  y = opt; 
  return
end
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
tic
xn=x;% remember the old file
[n m]= size(xn);
if n<m
 b=m;m=n;n=b; 
 xn=xn.';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

switch m
 case 1, xn=[(1:n)' xn(:)];
 case 2, % dimension OK.
 otherwise, 
   error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ')        
end

%%%%%%%%%%%%%%%%%%
%                %
%  initializing  %
%                %
%%%%%%%%%%%%%%%%%%
if nargin<2||isempty(inds),  inds=isnan(xn(:,2));end
if nargin<3||isempty(Nsim),  Nsim=20; end
if nargin<4||isempty(L),  L=[]; end, % lagsize
if nargin<5||isempty(def),  def='nonlinear';end
if nargin>=6,  opt=troptset(opt,varargin{:});end

param = opt.('param');

switch opt.plotflag
  case {'none','off'},   plotflag = 0;
  case 'final', plotflag = 1;
  case 'iter',  plotflag = 2;
  otherwise,    plotflag = opt.plotflag;
end

clf
olddef = def;
method = 'approx';  %'approx';%'dec2'; % 'dec2' 'approx';
ptime  = opt.delay; % pause for ptime sec if plotflag=2
pause on
expect1 = 1;     % first reconstruction by expectation? 1=yes 0=no
expect  = 1;     % reconstruct by expectation? 1=yes 0=no
tol     = 0.001; % absolute tolerance of e(g_new-g_old)

cmvmax = 100; % if number of consecutive missing values (cmv) are longer they
             % are not used in estimation of g, due to the fact that the
             % conditional expectation approaches zero as the length to
             % the closest known points increases, see below in the for loop 

dT=xn(2,1)-xn(1,1);
Lm=min([n,200,floor(200/dT)]); % Lagmax 200 seconds
if ~isempty(L), Lm=max([L,Lm]);end
Lma = 1500;                    % size of the moving average window used
                               % for detrending the reconstructed signal 



if any(inds>1),
  xn(inds,2)=NaN;
  inds=isnan(xn(:,2));
elseif sum(inds)==0,
  error('No spurious data given')
end
endpos  = diff(inds );
strtpos = find(endpos>0);
endpos  = find(endpos<0);

indg=find(~inds); % indices to good points
inds=find(inds);  % indices to spurous points


indNaN = []; % indices to points omitted in the covariance estimation
indr  = 1:n; % indices to point used in the estimation of g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finding more than cmvmax consecutive spurios points. 
% They will not be used in the estimation of g and are thus removed 
% from indr.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(strtpos) && (isempty(endpos)||endpos(end)<strtpos(end)),
  if (n-strtpos(end))>cmvmax 
    indNaN= indr(strtpos(end)+1:n);
    indr(strtpos(end)+1:n)=[];
  end
  strtpos(end)=[];
end
if ~isempty(endpos) && (isempty(strtpos)||(endpos(1)<strtpos(1))),
  if (endpos(1))>cmvmax 
    indNaN=[indNaN,indr(1:endpos(1))];
    indr(1:endpos(1))=[];
  end
  strtpos   = strtpos-endpos(1);
  endpos    = endpos-endpos(1);
  endpos(1) = [];
end
%length(endpos)
%length(strtpos)
for ix=length(strtpos):-1:1
  if (endpos(ix)-strtpos(ix)>cmvmax)
    indNaN=[indNaN, indr(strtpos(ix)+1:endpos(ix))];
    indr(strtpos(ix)+1:endpos(ix))=[]; % remove this when estimating the transform
  end
end
if length(indr)<0.1*n,
  error('Not possible to reconstruct signal')
end

if any(indNaN),
  indNaN=sort(indNaN);
end

switch 1, % initial reconstruction attempt
case  0,% spline
  xn(:,2) = interp1(xn(indg,1),xn(indg,2),xn(:,1),'*spline'); 
  y=xn;
  return
  
case 1,% 
 % xn(indg,2)=detrendma(xn(indg,2),1500);
  
  [g test cmax irr g2]  = dat2tr(xn(indg,:),def,opt);
  xnt=xn;
  xnt(indg,:)=dat2gaus(xn(indg,:),g);
  xnt(inds,2)=NaN;
  rwin=findrwin(xnt,Lm,L);
  disp(['First reconstruction attempt,    e(g-u)=', num2str(test)] )
  [samp ,mu1o, mu1oStd]  = cov2csdat(xnt(:,2),rwin,1,method,inds); % old simcgauss
  if expect1,% reconstruction by expectation
    xnt(inds,2) =mu1o;
  else
    xnt(inds,2) =samp;
  end
  xn=gaus2dat(xnt,g);
  xn(:,2)=detrendma(xn(:,2),Lma); % detrends the signal with a moving
                                   % average of size Lma
  g_old=g;
  
end

bias = mean(xn(:,2));
xn(:,2)=xn(:,2)-bias; % bias correction

if plotflag==2
  clf
  mind=1:min(1500,n);
  waveplot(xn(mind,:),x(inds(mind),:), 6,1)
  subplot(111)
  pause(ptime)
end

test0=0;
for ix=1:Nsim,
%   if 0,%ix==2,
%     rwin=findrwin(xn,Lm,L);
%     xs=cov2sdat(rwin,[n 100 dT]);
%     [g0 test0 cmax irr g2]  = dat2tr(xs,def,opt);
%     [test0 ind0]=sort(test0);
%   end
  
   if 1, %test>test0(end-5),
     % 95% sure the data comes from a non-Gaussian process
     def = olddef; %Non Gaussian process
   else
     def = 'linear'; % Gaussian process
   end
   % used for isope article
   % indr =[1:27000 30000:39000];
   % Too many consecutive missing values will influence the estimation of
   % g. By default do not use consecutive missing values if there are more 
   % than cmvmax. 
   
   [g test cmax irr g2]  = dat2tr(xn(indr,:),def,opt);
  if plotflag==2,
    pause(ptime)
  end
  
  
  %tobs=sqrt((param(2)-param(1))/(param(3)-1)*sum((g_old(:,2)-g(:,2)).^2))
  % new call
  tobs=sqrt((param(2)-param(1))/(param(3)-1)....
	*sum((g(:,2)-interp1(g_old(:,1)-bias, g_old(:,2),g(:,1),'spline')).^2));
  
  if ix>1 
    if tol>tobs2 && tol>tobs,	
      break, %estimation of g converged break out of for loop    
    end
  end
 
  tobs2=tobs;
  
  xnt=dat2gaus(xn,g);
  if ~isempty(indNaN),    xnt(indNaN,2)=NaN;  end
  rwin=findrwin(xnt,Lm,L);    
  disp(['Simulation nr: ', int2str(ix), ' of ' num2str(Nsim),'   e(g-g_old)=', num2str(tobs), ',  e(g-u)=', num2str(test)])
  [samp ,mu1o, mu1oStd]  = cov2csdat(xnt(:,2),rwin,1,method,inds);
  
  if expect,
    xnt(inds,2) =mu1o;
  else
    xnt(inds,2) =samp;
  end
  
  xn=gaus2dat(xnt,g);
  if ix<Nsim
    bias=mean(xn(:,2));
    xn(:,2) = (xn(:,2)-bias); % bias correction
  end
  g_old=g;% saving the last transform
  if plotflag==2
    waveplot(xn(mind,:),x(inds(mind),:),6,1,[])
    subplot(111)
    pause(ptime)
  end
end % for loop

if 1, %test>test0(end-5) 
  xnt=dat2gaus(xn,g);
  [samp ,mu1o, mu1oStd]  = cov2csdat(xnt(:,2),rwin,1,method,inds);
  xnt(inds,2) =samp;
  xn=gaus2dat(xnt,g);
  bias=mean(xn(:,2));
  xn(:,2) = (xn(:,2)-bias); % bias correction
  g(:,1)=g(:,1)-bias;
  g2(:,1)=g2(:,1)-bias;
  gn=trangood(g);
 
  %mu1o=mu1o-tranproc(bias,gn);
  muUStd=tranproc(mu1o+2*mu1oStd,fliplr(gn));%
  muLStd=tranproc(mu1o-2*mu1oStd,fliplr(gn));%
else
  muLStd=mu1o-2*mu1oStd;
  muUStd=mu1o+2*mu1oStd;
end

if  plotflag==2 && length(xn)<10000,
  waveplot(xn,[xn(inds,1) muLStd ;xn(inds,1) muUStd ], 6,round(n/3000),[])
  legend('reconstructed','2 stdev')
  %axis([770 850 -1 1])
  %axis([1300 1325 -1 1])
end
y=xn;
toc

return

function r=findrwin(xnt,Lm,L)
  r=dat2cov(xnt,Lm);%computes  ACF
  %finding where ACF is less than 2 st. deviations .
  % in order to find a better L  value
  if nargin<3||isempty(L)
    L=find(abs(r.R)>2*r.stdev)+1;
    if isempty(L), % pab added this check 09.10.2000
      L = Lm;
    else
      L = min([floor(4/3*L(end)) Lm]);
    end
  end
  win=parzen(2*L-1);
  r.R(1:L)=win(L:2*L-1).*r.R(1:L);
  r.R(L+1:end)=0;
  return
