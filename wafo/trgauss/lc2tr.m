function [g, test, g2] = lc2tr(cross,ma,sa,varargin)
%LC2TR Estimate transformation, g, from observed crossing intensity.
%
%        Assumption: a Gaussian process, Y, is related to the
%                    non-Gaussian process, X, by Y = g(X). 
%       
%  CALL:  [g, test,g2] = lc2tr(lc,ma,sa,options);
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
% csm, gsm  - defines the smoothing of the logarithm of crossing intensity 
%             and the transformation g, respectively. Valid values must 
%             be 0<=csm,gsm<=1. (default csm=0.9, gsm=0.05)
%             Smaller values gives smoother functions.
%  param    - vector which defines the region of variation of the data X.
%             (default [-5 5 513]). 
%  plotflag - 0 no plotting (Default)
%             1 plots empirical and smoothed g(u) and the theoretical for a Gaussian model. 
%             2 monitor the development of the estimation
% linextrap - 0 use a regular smoothing spline 
%             1 use a smoothing spline with a constraint on the ends to 
%               ensure linear extrapolation outside the range of the data.
%               (default)
% cvar      - Variances for the logarithm of the crossing intensity. (default  1) 
% gvar      - Variances for the empirical transformation, g. (default  1) 
% ne        - Number of extremes (maxima & minima) to remove from the
%              estimation of the transformation. This makes the
%              estimation more robust against outliers. (default 7)
%
%    The empirical crossing intensity is usually very irregular.
%  More than one local maximum of the smoothed crossing intensity
%  may cause poor fit of the transformation. In such case one
%  should use a smaller value of CSM or set a larger variance for CVAR. 
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
% g0 = lc2tr(lc,0,Hm0/4,'plot','iter');         % Monitor the development
% g1 = lc2tr(lc,0,Hm0/4,troptset('gvar', .5 )); % Equal weight on all points
% g2 = lc2tr(lc,0,Hm0/4,'gvar', [3.5 .5 3.5]);  % Less weight on the ends
% hold on, trplot(g1,g)                         % Check the fit
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
% revised pab 21.12.2000
% - added interp of cross -> much faster estimation
% - added example, chkder
% - replaced optional arguments with a options struct
% - added options to input: ne, cvar,gvar
% - changed names: csm1 to csm
%                  csm2 to gsm
% - replaced monitor with plotflag==2
% - default param is now [-5 5 513] -> better to have the discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.
%
%  modified by svi 29.09.99
%  Transformations  g2, g2 are not normalized any longer.
% revised pab 11.08.99
% changed name from cross2tr to lc2tr
%
% modified by Per A. Brodtkorb 15.08.98
%  to check if the smoothing is sufficient and 
%  changed the calculation of the test statistic.
%  moved the plotting routine to trplot 

%  Beware of the problem that Carl de Boor's smooth function 
%  does not always extrapolate well outside the ends when 
%  the smoothing parameter, p, is close to one.
%  Particularly extrapolation in the first smoothing may corrupt 
%  the estimate of g. One solution to the problem is to 
%  extrapolate linearly. This is incorporated into the smooth function.
%  Yet, another solution is choosing a lower value for csm1
%  or not to extrapolate at all in the first smoothing but instead leave
%  all the extrapolation to the second smoothing. 
%  (Probably better since csm2<<csm1))
%
%   also added secret options: plotflag and monitor the steps of
%   estimation of the transformation 
%

opt = troptset('chkder','on','plotflag','off','csm',.95,'gsm',.05,....
    'param',[-5 5 513],'delay',2,'ntr',1000,'linextrap','on','ne',7,'cvar',1,'gvar',1);
% If just 'defaults' passed in, return the default options in g
if nargin==1 && nargout <= 1 && isequal(cross,'defaults')
  g = opt; 
  return
end
%error(nargchk(3,inf,nargin)) 
narginchk(3,inf)
if nargin>=4,  opt  = troptset(opt,varargin{:}); end

csm1 = opt.csm;
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
   %cros = [ x0,interp1q(cross(:,1),cross(:,2),x0)];
   Ne=0;
   ncr = opt.ntr;
else
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
cros(:,1) = (cros(:,1)-ma)/sa;
%cros0=cros;

if 0, % slightly smoothing the crossing spectrum
  tmp=cssmooth(cros(1:end,1),cros(1:end,2),0.95,cros(1:end,1)); 
  ind=(tmp>0);
  cros(ind,2)=tmp; clear tmp ind
end

indz = (cros(:,2)==0);
if any(indz)
  cros(indz,2) = inf; % this is done in order to avoid a warning message
end
cros(~indz,2) = -log(cros(~indz,2));


% NB! the cssmooth function does not always extrapolate well outside the edges
% causing poor estimate of g  
% We may alleviate this problem by: forcing the extrapolation
% to be linear outside the edges or choosing a lower value for csm1
% or not to extrapolate at all in the first smoothing but instead
% extrapolate in the second smoothing. (Possibly better since csm2<<csm1)
% Therefore replacing the old call
%scros=cssmooth(cros(10:ncr-10,1),cros(10:ncr-10,2),csm1,cros(1:ncr,1)); 
% with
inds = 1+Ne:ncr-Ne;% indices to points we are smoothing over
inde = 1+Ne:ncr-Ne;% indices to points we are smoothing over and
             % possibly  extrapolating if length(inds)<length(inde)
scros = cssmooth(cros(inds,1),cros(inds,2),csm1,cros(inde,1),def,cvar(inds)); 

if plotflag>1
  %plot(cros(10:ncr-10,1),cros(10:ncr-10,2),'r')
  plot(cros(:,1),cros(:,2),'r'),
  hold on
  plot(cros(inde,1),scros,'b'),hold off
  title('First smoothing')
  ylabel('-log(cross)')
  xlabel('crossing level')
  legend('Not smoothed','Smoothed',0)
  %return
  pause(ptime)       
end

%  scros has to have a single minimum
if 0,% old call
  [smin imin]=min(scros);%
else % new call: checking if we have a single minimum 
 imin = findcross(diff(scros))+1; 
 smin = scros(imin); 
 if length(imin)~=1,
   disp(['Warning:  There are ' num2str(length(imin)) ' minima/' ...
	   'maxima after the first smoothing'])  
   [smin ind] = min(smin);
   imin       = imin(ind);
 end
end
%imin=imin+30
%smin = scros(imin)

scros1 = sqrt(2*abs(scros-smin));
%scros1(1:imin)=-scros1(1:imin);
scros1(1:imin) = 2*scros1(imin)-scros1(1:imin);

scros2 = cssmooth(cros(inde,1),scros1,csm2,uu,def,gvar(inde));

g(:,2) = scros2';%*sa; %multiply with stdev 

if nargout>2||plotflag>0,
  cros2 = cros; 
  [cmin icmin] = min(cros(:,2));
  cros2(:,2)   = sqrt(2*abs(cros(:,2)-cmin));
  cros2(1:icmin,2) = 2*cros2(icmin,2)-cros2(1:icmin,2); 
  g2(:,2)  = cros2(:,2);
end

if chkder~=0
   for ix = 1:5
    dy = diff(g(:,2));
    if any(dy<=0)
      warning('WAFO:LC2TR','The empirical crossing spectrum is not sufficiently smoothed.')
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
% wrong!!!                
%test=0.02*sqrt(sum((uu-scros2).^2));
%                  5 
%this is not sqrt(int (g(u)-u)^2 du)
%                 -5
% The correct is
if 0,  %either
  test = sqrt((param(2)-param(1))/(param(3)-1)*sum((uu-scros2).^2));
else % or
  %test=sqrt(simpson(uu,(uu-scros2).^2));
% or
  test=sqrt(trapz(uu,(uu-scros2).^2));
end 


if plotflag>0,
  %% Plotchanges made by Joakim Elvander 970707.
  %% Plots will be initiated in the file funplot_1
  
  % %clf
  % plot(uu,scros2,'r')
  % hold on
  % plot(uu,uu,'g--')
  % pause
  
  
  if plotflag>1
    stairs(cros2(:,1),cros2(:,2),'b'),hold on
    plot(cros(inde,1),scros1,'r')
    plot(uu,scros2,'g--'),hold off
    title('Second smoothing')
    ylabel('Transfer function g(u)')
    xlabel('Crossing level u')
    legend('Not smoothed','Smoothed once','Smoothed twice',0)
    pause(ptime)
  end
  
  
  %stairs(cros2(:,1),cros2(:,2))  
  %axis([-5 5 -5 5])
  %axis('square')
  %hold off
 
  %% Temporarily
  trplot(g,g2,ma,sa)
  %  funplot_1(uu,scros2,cros2);
  %legend(['Smoothed (T='  num2str(test) ')'],'g(u)=u','Not smoothed',0)
  %ylabel('Transfer function g(u)')
  %xlabel('Crossing level u')
  
  if plotflag>1,pause(ptime),end
  %hold on, stairs(cros2(10:ncr-10,1),cros2(10:ncr-10,2),'y');hold off
end







