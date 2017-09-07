function [g, test, g2] = lc2tr2(cross,ma,sa,varargin)
%LC2TR2 Estimate transformation, g, from observed crossing intensity, version2.
%
%        Assumption: a Gaussian process, Y, is related to the
%                    non-Gaussian process, X, by Y = g(X). 
%
%  CALL:  [g, test,g2] = lc2tr2(lc,ma,sa,options);
%
%     g,g2  = smoothed and empirical estimate of the transformation  g.     
%     test  = test observator int (g(u)-u)^2 du  where int limits is
%             given by param. This is a measure of departure of the 
%             data from the Gaussian model.
%     lc    = a two column matrix with levels in the first column
%             and number of upcrossings in the second.
%     ma,sa = mean and standard deviation of the process
%
%   options = structure with the fields:
%  csm,gsm  - defines the smoothing of the crossing intensity and the
%             transformation g. Valid values must be
%             0<=csm,gsm<=1. (default csm = 0.9 gsm=0.05)
%             Smaller values gives smoother functions.
%     param - vector which defines the region of variation of the data X.
%             (default [-5 5 513]). 
%  plotflag - 0 no plotting (Default)
%             1 plots empirical and smoothed g(u) and theoretical for a Gaussian model. 
%             2 monitor development of estimation
% linextrap - 0 use a regular smoothing spline 
%             1 use a smoothing spline with a constraint on the ends to 
%               ensure linear extrapolation outside the range of the data.
%               (default)
% cvar      - Variances for the crossing intensity. (default 1)
% gvar      - Variances for the empirical transformation, g. (default  1) 
% ne        - Number of extremes (maxima & minima) to remove from the
%              estimation of the transformation. This makes the
%              estimation more robust against outliers. (default 7)
% Ntr        - Maximum length of empirical crossing intensity.
%              The empirical crossing intensity is interpolated
%              linearly  before smoothing if the length exceeds Ntr.
%              A reasonable NTR will significantly speed up the
%              estimation for long time series without loosing any
%              accuracy. NTR should be chosen greater than
%              PARAM(3). (default 1000)
%
%    The empirical crossing intensity is usually very irregular.
%  More than one local maximum of the empirical crossing intensity
%  may cause poor fit of the transformation. In such case one
%  should use a smaller value of GSM or set a larger variance for GVAR. 
%    If X(t) is likely to cross levels higher than 5 standard deviations  
%  then the vector param has to be modified.  For example if X(t) is 
%  unlikely to cross a level of 7 standard deviations one can use 
%  param = [-7 7 513].
%
% Example:
% Hm0 = 7;
% S = jonswap([],Hm0); g=ochitr([],[Hm0/4]); 
% S.tr=g;S.tr(:,2)=g(:,2)*Hm0/4;
% xs = spec2sdat(S,2^13);
% lc = dat2lc(xs);
% g0 = lc2tr2(lc,0,Hm0/4,'plot','iter');         % Monitor the development
% g1 = lc2tr2(lc,0,Hm0/4,troptset('gvar', .5 )); % Equal weight on all points
% g2 = lc2tr2(lc,0,Hm0/4,'gvar', [3.5 .5 3.5]);  % Less weight on the ends
% hold on, trplot(g1,g)                          % Check the fit
% trplot(g2)
%
% See also  troptset, dat2tr, trplot, findcross, cssmooth

% NB! the transformated data will be N(0,1)

% Reference
% Rychlik , I., Johannesson, P., and Leadbetter, M.R. (1997)
% "Modelling and statistical analysis of ocean wavedata 
% using a transformed Gaussian process",
% Marine structures, Design, Construction and Safety, 
% Vol 10, pp 13--47
% 

% Adapted to  cssmooth  by GL Feb 2011
% Tested on: Matlab 5.3, 5.2, 5.1
% History:
% by pab 29.12.2000
% based on lc2tr, but the inversion is faster.
% by IR and PJ


opt = troptset('chkder','on','plotflag','off','csm',.9,'gsm',.05,....
    'param',[-5 5 513],'delay',2,'linextrap','on','ntr',1000,'ne',7,'cvar',1,'gvar',1);
% If just 'defaults' passed in, return the default options in g
if nargin==1 && nargout <= 1 && isequal(cross,'defaults')
  g = opt; 
  return
end
%error(nargchk(3,inf,nargin)) 
narginchk(3,inf)
if nargin>=4 ,  opt  = troptset(opt,varargin{:}); end
csm2 = opt.gsm;
param = opt.param;
ptime = opt.delay;
Ne  = opt.ne;
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
ncr = length(cross);
if ncr>opt.ntr && opt.ntr>0,
   x0 = linspace(cross(1+Ne,1),cross(end-Ne,1),opt.ntr)';
   cros = [ x0,interp1(cross(:,1),cross(:,2),x0, 'linear')];
   % cros = [ x0,interp1q(cross(:,1),cross(:,2),x0)];
   Ne = 0;
   Ner = opt.ne;
   ncr = opt.ntr;
 else
   Ner = 0;
  cros=cross;
end

ng = length(opt.gvar);
if ng==1
  gvar = opt.gvar(ones(ncr,1));
else
  gvar = interp1(linspace(0,1,ng)',opt.gvar(:),linspace(0,1,ncr)','*linear');  
end
ng = length(opt.cvar);
if ng==1
  cvar = opt.cvar(ones(ncr,1));
else
  cvar = interp1(linspace(0,1,ng)',opt.cvar(:),linspace(0,1,ncr)','*linear');  
end

g = zeros(param(3),2);

uu = levels(param);

g(:,1) = sa*uu' + ma;

g2   = cros;

if Ner>0, % Compute correction factors
 cor1 = trapz(cross(1:Ner+1,1),cross(1:Ner+1,2));
 cor2 = trapz(cross(end-Ner-1:end,1),cross(end-Ner-1:end,2));
else
  cor1 = 0;
  cor2 = 0;
end
cros(:,2) = cumtrapz(cros(:,1),cros(:,2))+cor1;
cros(:,2) = (cros(:,2)+.5)/(cros(end,2) + cor2 +1);
cros(:,1) = (cros(:,1)-ma)/sa;

% find the mode
[tmp,imin]= min(abs(cros(:,2)-.15));
[tmp,imax]= min(abs(cros(:,2)-.85));
inde = imin:imax;
tmp =  cssmooth(cros(inde,1),g2(inde,2),opt.csm,cros(inde,1),def,cvar(inde));

[tmp imax] = max(tmp);
u0 = cros(inde(imax),1);
%u0 = interp1q(cros(:,2),cros(:,1),.5)


cros(:,2) = invnorm(cros(:,2),-u0,1);

g2(:,2)   = cros(:,2);
% NB! the smooth function does not always extrapolate well outside the edges
% causing poor estimate of g  
% We may alleviate this problem by: forcing the extrapolation
% to be linear outside the edges or choosing a lower value for csm2.

inds = 1+Ne:ncr-Ne;% indices to points we are smoothing over
scros2 = cssmooth(cros(inds,1),cros(inds,2),csm2,uu,def,gvar(inds));

g(:,2) = scros2';%*sa; %multiply with stdev 

if chkder~=0
   for ix = 1:5
    dy = diff(g(:,2));
    if any(dy<=0)
      warning('WAFO:LCTR2','The empirical crossing spectrum is not sufficiently smoothed.')
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
if 0,  %either
  test = sqrt((param(2)-param(1))/(param(3)-1)*sum((uu-scros2).^2));
else % or
  %test=sqrt(simpson(uu,(uu-scros2).^2));
% or
  test=sqrt(trapz(uu,(uu-scros2).^2));
end 


if plotflag>0, 
  trplot(g ,g2,ma,sa)
  %legend(['Smoothed (T='  num2str(test) ')'],'g(u)=u','Not smoothed',0)
  %ylabel('Transfer function g(u)')
  %xlabel('Crossing level u')
  
  if plotflag>1,pause(ptime),end
end







