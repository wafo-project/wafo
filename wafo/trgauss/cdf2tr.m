function [g, test, g2] = cdf2tr(Fx1,ma ,sa,varargin)
%CDF2TR Estimate transformation, g, from observed CDF.
%
%        Assumption: a Gaussian process, Y, is related to the
%                    non-Gaussian process, X, by Y = g(X). 
% 
%  CALL [g,test,g2] = cdf2tr(F,ma,sa,options);
%
%     g,g2  = smoothed and empirical estimate of the transformation  g.     
%     test  = test observator int (g(u)-u)^2 du  where int limits is
%             given by OPTIONS.PARAM. This is a measure of departure of the 
%             data from the Gaussian model.
%     F     = empirical CDF of X(t), a 2 column matrix.
%     ma,sa = mean and standard deviation of the process X(t).
%   options = options structure defining how the smoothing is done.
%             (See troptset for default values)
%
%    The empirical CDF is usually very irregular.
%  More than one local maximum of the empirical CDF
%  may cause poor fit of the transformation. In such case one
%  should use a smaller value of GSM or set a larger variance for GVAR. 
%    If X(t) is likely to cross levels higher than 5 standard deviations  
%  then the vector param has to be modified.  For example if X(t) is 
%  unlikely to cross a level of 7 standard deviations one can use 
%  param = [-7 7 513].
%
% Example
% Hm0 = 7;
% S = jonswap([],Hm0); g=ochitr([],[Hm0/4]); 
% S.tr = g; S.tr(:,2)=g(:,2)*Hm0/4;
% xs = spec2sdat(S,2^13);
% Fx = edf(xs(:,2),'wdata',false);
% lc = dat2lc(xs);
% g0 = cdf2tr(Fx,0,Hm0/4,troptset('plot',1));   % Plot final estimate
% g1 = lc2tr(lc,0,Hm0/4,troptset('gvar', .5 )); % More weight on all points
% g2 = lc2tr(lc,0,Hm0/4,'gvar', [3.5 .5 3.5]);  % Less weight on the ends
% hold on, trplot(g1,g)                         % Check the fit
% trplot(g2)
%
% See also  troptset, lc2tr

% Adapted to  cssmooth  by GL Feb 2011
% History
% revised Feb2004  
% Revised pab Dec2003
%  fixed a bug F -> Fx1
% by pab 29.12.2000
% - default param is now [-5 5 513] -> better to have the discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.

opt = troptset('chkder','on','plotflag','off','gsm',.05,....
    'param',[-5 5 513],'delay',2,'linextrap','on','ntr',1000,'ne',7,'gvar',1);
% If just 'defaults' passed in, return the default options in g
if nargin==1 && nargout <= 1 && isequal(Fx1,'defaults')
  g = opt; 
  return
end
%error(nargchk(3,inf,nargin))
narginchk(3,inf)
if nargin>=4,  opt=troptset(opt,varargin{:}); end
switch opt.chkder;
  case 'off', chkder = 0;
  case 'on',  chkder = 1;
  otherwise,  chkder = opt.chkder;
end
switch opt.linextrap;
  case 'off', def = 0;
  case 'on',  def = 1;
  otherwise,  def = opt.linextrap;
end

switch opt.plotflag
  case {'none','off'},   plotflag = 0;
  case 'final', plotflag = 1;
  case 'iter',  plotflag = 2;
  otherwise,    plotflag = opt.plotflag;
end


Ne = opt.ne;
if length(Fx1)>opt.ntr && opt.ntr>0
  x0 = linspace(Fx1(1+Ne,1),Fx1(end-Ne,1),opt.ntr)';
  Fx = [ x0,interp1(Fx1(:,1),Fx1(:,2),x0, 'linear')];
  % Fx = [ x0,interp1q(Fx1(:,1),Fx1(:,2),x0)];
  Ne=0;
else
  Fx = Fx1;
end
uu = levels(opt.param)';
g = [sa*uu+ma zeros(opt.param(3),1)];
ncr = length(Fx);


ng = length(opt.gvar);
if ng==1
  gvar = opt.gvar(ones(ncr,1));
else
  gvar = interp1(linspace(0,1,ng)',opt.gvar(:),linspace(0,1,ncr)','*linear');  
end
 

ind = find(diff(Fx(:,1))>0); % remove equal points
ind1 = ind(Ne+1:end-Ne);  
tmp = invnorm(Fx(ind,2));

g(:,2) = cssmooth(Fx(ind1,1),tmp(Ne+1:end-Ne),opt.gsm,g(:,1),def,gvar);

if chkder~=0
  for ix = 1:5
    dy = diff(g(:,2));
    if any(dy<=0)
      warning('WAFO:CDF2TR','The empirical distribution is not sufficiently smoothed.')
      disp('        The estimated transfer function, g, is not ')
      disp('        a strictly increasing function.')
      dy(dy>0)=eps;
      gvar = -([dy;0]+[0;dy])/2+eps;
      g(:,2) = cssmooth(g(:,1),g(:,2),1,g(:,1),def,ix*gvar);
    else 
      break
    end
  end
end
if nargout>1
  test = sqrt(trapz(uu,(uu-g(:,2)).^2));  
end
if nargout>2 || plotflag>0
  g2 = [Fx(ind,1) tmp];
end
if plotflag>0
  trplot(g,g2,ma,sa)
  if plotflag>1
    pause(opt.delay)
  end
end



