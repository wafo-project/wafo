function [g, test, cmax, irr, g2]= dat2tr(x,def,varargin)
%DAT2TR Estimate transformation, g, from data.
%
% CALL:  [g test cmax irr g2]  = dat2tr(x,def,options);
%
%   g,g2   = the smoothed and empirical transformation, respectively. 
%            A two column matrix if multip=0.  
%            If multip=1 it ís a 2*(m-1) column matrix where the
%            first and second column is the transform 
%            for values in column 2 and third and fourth column is the
%            transform for values in column 3 ......
%
%   test   = int (g(u)-u)^2 du  where int. limits is given by param. This
%            is a measure of departure of the data from the Gaussian model.
%           
%   cmax   = maximum crossing intensity of x
%   irr    = irregularity factor of x which is approximately Tz/Tmaxima   
%   x      = m column data matrix with sampled times in the first column
%            and values the next columns.            
%
%   def    = 'nonlinear' : transform based on smoothed crossing intensity (default)
%            'mnonlinear': transform based on smoothed marginal distribution
%            'hermite'   : transform based on cubic Hermite polynomial
%            'ochitr'    : transform based on exponential function
%            'linear'    : identity.
%
%   options = options structure with the following fields:
%  csm,gsm - defines the smoothing of the logarithm of crossing intensity 
%            and the transformation g, respectively. Valid values must 
%            be 0<=csm,gsm<=1. (default csm=0.9, gsm=0.05)
%            Smaller values gives smoother functions.
%    param - vector which defines the region of variation of the data x.
%           (default see lc2tr). 
% plotflag - 0 no plotting (Default)
%            1 plots empirical and smoothed g(u) and the theoretical for
%              a Gaussian model. 
%            2 monitor the development of the estimation
%linextrap - 0 use a regular smoothing spline 
%            1 use a smoothing spline with a constraint on the ends to 
%              ensure linear extrapolation outside the range of the data.
%              (default)
%     gvar - Variances for the empirical transformation, g. (default  1) 
%       ne - Number of extremes (maxima & minima) to remove from the
%            estimation of the transformation. This makes the
%            estimation more robust against outliers. (default 7)
%      ntr - Maximum length of empirical crossing intensity or CDF.
%            The empirical crossing intensity or CDF is interpolated
%            linearly  before smoothing if their lengths exceeds Ntr.
%            A reasonable NTR will significantly speed up the
%            estimation for long time series without loosing any
%            accuracy. NTR should be chosen greater than
%            PARAM(3). (default 1000)
%   multip - 0 the data in columns belong to the same seastate (default).
%            1 the data in columns are from separate seastates.
%
%  DAT2TR estimates the transformation in a transformed Gaussian model.  
%  Assumption: a Gaussian process, Y, is related to the
%  non-Gaussian process, X, by Y = g(X). 
% 
%  The empirical crossing intensity is usually very irregular.
%  More than one local maximum of the empirical crossing intensity
%  may cause poor fit of the transformation. In such case one
%  should use a smaller value of CSM. In order to check the effect 
%  of smoothing it is recomended to also plot g and g2 in the same plot or
%  plot the smoothed g against an interpolated version of g (when CSM=GSM=1).
%    If  x  is likely to cross levels higher than 5 standard deviations
%  then the vector param has to be modified.  For example if x is 
%  unlikely to cross a level of 7 standard deviations one can use 
%  PARAM=[-7 7 513].
%
% Example:
% Hm0 = 7;
% S = jonswap([],Hm0); g=ochitr([],[Hm0/4]); 
% S.tr=g;S.tr(:,2)=g(:,2)*Hm0/4;
% xs = spec2sdat(S,2^13);
% g0 = dat2tr(xs,[],'plot','iter');             % Monitor the development
% g1 = dat2tr(xs,'mnon','gvar', .5 );           % More weight on all points
% g2 = dat2tr(xs,'nonl','gvar', [3.5 .5 3.5]);  % Less weight on the ends
% hold on, trplot(g1,g)                                   % Check the fit
% trplot(g2)
%
% See also  troptset, lc2tr, cdf2tr, trplot

% References:
% Rychlik, I. , Johannesson, P and Leadbetter, M. R. (1997)
% "Modelling and statistical analysis of ocean wavedata using 
%  transformed Gaussian process."
% Marine structures, Design, Construction and Safety, Vol. 10, No. 1, pp 13--47
%
% 
% Brodtkorb, P, Myrhaug, D, and Rue, H (1999)
% "Joint distribution of wave height and crest velocity from
% reconstructed data"
% in Proceedings of 9th ISOPE Conference, Vol III, pp 66-73



%Tested on: Matlab 5.3, 5.2, 5.1
%History:
% revised pab Dec2004
%  -Fixed bug: string comparison for def at fault.  
% revised pab Nov2004
%  -Fixed bug: linextrap was not accounted for  
% revised pab july 2004
% revised pab 3 april 2004
% -fixed a bug in hermite estimation: excess changed to kurtosis  
% revised pab 29.12.2000
% - added example, hermite and ochi options
% - replaced optional arguments with a options struct
% - default param is now [-5 5 513] -> better to have the discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.
% revised pab 19.12.2000
%  - updated call edf(X,-inf,[],monitor) to  edf(X,[],monitor)
%    due to new calling syntax for edf
% modifed pab 24.09.2000
%  - changed call from norminv to wnorminv
%  - also removed the 7 lowest and 7 highest points from
%    the estimation using def='mnonlinear' 
%    (This is similar to what lc2tr does. lc2tr removes
%     the 9 highest and 9 lowest TP from the estimation)
% modified pab 09.06.2000
%  - made all the *empirical options secret.
%  - Added 'mnonlinear' and 'mempirical' 
%  - Fixed the problem of multip==1 and def=='empirical' by interpolating 
%    with spline to ensure that the length of g is fixed
%  - Replaced the test statistic for def=='empirical' with the one
%    obtained when csm1=csm2=1. (Previously only the smoothed test
%    statistic where returned)
% modified pab 12.10.1999
%  fixed a bug
%  added secret output of empirical estimate g2
% modified by svi  29.09.1999
% changed input def by adding new options.
% revised by pab 11.08.99
%   changed name from dat2tran to dat2tr
% modified by Per A. Brodtkorb 12.05.1999,15.08.98
%   added  secret option: to accept multiple data, to monitor the steps 
%   of estimation of the transformation 
%   also removed some code and replaced it with a call to lc2tr (cross2tr) 
%   making the maintainance easier
%

%opt = troptset('plotflag','off','csm',.95,'gsm',.05,....
%    'param',[-5 5 513],'delay',2,'linextrap','on','ne',7,...
%    'cvar',1,'gvar',1,'multip',0);


opt = troptset('chkder','on','plotflag','off','csm',.95,'gsm',.05,....
    'param',[-5 5 513],'delay',2,'ntr',1000,'linextrap','on','ne',7,'cvar',1,'gvar',1,'multip',0,'crossdef','uM');
% If just 'defaults' passed in, return the default options in g
if nargin==1 && nargout <= 1 && isequal(x,'defaults')
  g = opt; 
  return
end
error(nargchk(1,inf,nargin)) 
if nargin<2||isempty(def),     def    = 'nonlinear'; end
if nargin>=3,  opt   = troptset(opt,varargin{:});   end
multip = opt.multip;
%Ne = opt.ne;
switch opt.plotflag
  case {'none','off'},   plotflag = 0;
  case 'final', plotflag = 1;
 case 'iter',  plotflag = 2;
  
  otherwise,    plotflag = opt.plotflag;
end

% Crossing definition
switch opt.crossdef
  case {'u'}, cdef = 1; % only upcrossings.
  case 'uM',  cdef = 2; % upcrossings and Maxima (default).
  case 'umM', cdef = 3; % upcrossings, minima, and Maxima.
  case 'um',  cdef = 4; % upcrossings and minima.  
  otherwise,  cdef = opt.crossdef;
end
switch opt.linextrap
 case {'on'},  linextrap = 1;
 case {'off'}, linextrap = 0;
 otherwise,    linextrap = opt.linextrap;
end

xx = x;
[n,m] = size(xx);
ma = mean(xx(:,2:m));
sa = std(xx(:,2:m)); 
m2 = m;

if m>2 && multip==0,% data in columns belongs to the same seastate
  ma = mean(ma);
  sa = sqrt(mean(sa.^2)); % pooled standard deviation
  m2 = 2;
end

  
g  = zeros(opt.param(3),2*(m2-1));
uu = levels(opt.param)';
g(:,2:2:end) = uu(:,ones(1,m2-1));
g(:,1:2:end) = sa(ones(opt.param(3),1),:).*g(:,2:2:end) + ...
      ma(ones(opt.param(3),1),:);
g2 = g;

test = zeros(m2-1,1);

if strncmpi(def,'lin',3) && nargout<=2, return,end

irr   = test;
cmax  = irr;


if multip==1,
  for ix=1:(m-1),
    if (lower(def(1))=='n'||lower(def(1))=='e' || nargout>2),
      tp       = dat2tp(xx(:,[1 ix+1]));% Turning points 
      mM       = tp2mm(tp);             % min2Max cycles
      cross1   = mm2lc(mM,cdef,0);      % Want upcrossings and maxima
      cmax(ix) = max(cross1(:,2));
      irr(ix)  = length(mM)/cmax(ix);   % approximately Tz/Tmaxima
    end
    if (lower(def(1))=='m')
      Fx  = edf(xx(:, ix+1),'wdata',false);
      if plotflag==2
        stairs(Fx(:,1),Fx(:,2))
        pause(opt.delay)
      end
    end
    switch lower(def(1:3))
      case {'her'},
        ga1 = skew(xx(:,ix+1));
        ga2 = kurt(xx(:,ix+1))-3;
        up  = min(4*(4*ga1/3).^2,12);
        lo  = ga1^2*3/2;
        kurt1 = min(up,max(ga2,lo))+3;
        phat = [sa(ix), ga1, kurt1,ma(ix) ];
        [g(:,2*ix-1:2*ix), test(ix)] = hermitetr(g(:,2*ix-1) ,phat);
        g2=[];
      case {'och'},
        phat = [ sa(ix) skew(xx(:,ix+1)) ma(ix)];
        [g(:,2*ix-1:2*ix), test(ix)] = ochitr(g(:,2*ix-1) ,phat);
        g2=[];
      case {'non'}, % nonlinear	  
        [g(:,2*ix-1:2*ix), test(ix),tmp]=lc2tr(cross1,ma(ix),sa(ix),opt);
        if nargout>4
          g2(:,2*ix-1) = g(:,2*ix-1);
          g2(:,2*ix)   = smooth(tmp(:,1),tmp(:,2),1,g(:,2*ix-1),linextrap);
        end
      case {'emp'}, % empirical 
        [g2(:,2*ix-1:2*ix), test(ix),tmp]=lc2tr(cross1,ma(ix),sa(ix),opt);
	g(:,2*ix-1) = g2(:,2*ix-1);
	g(:,2*ix)   = smooth(tmp(:,1),tmp(:,2),1,g2(:,2*ix-1),1);
	test(ix)    = sqrt(trapz(uu,(uu-g(:,2*ix)).^2));
      case {'mno'} % mnonlinear 
	[g(:,2*ix-1:2*ix), test(ix),tmp]= cdf2tr(Fx,ma(ix),sa(ix),opt);
	if nargout>4
	  g2(:,2*ix-1) = g(:,2*ix-1);
	  g2(:,2*ix)   = smooth(tmp(:,1),tmp(:,2),1,g(:,2*ix-1),1);
	end
      case {'mem'}, % mempirical
	[g2(:,2*ix-1:2*ix), test(ix),tmp]=cdf2tr(Fx,ma(ix),sa(ix),opt);
	g(:,2*ix-1) = g2(:,2*ix-1);
	g(:,2*ix)   = smooth(tmp(:,1),tmp(:,2),1,g2(:,2*ix-1),linextrap);
	test(ix)    = sqrt(trapz(uu,(uu-g(:,2*ix)).^2));
	%g(:,2*ix)  = smooth(Fx(ind,1),tmp,1,g(:,2*ix-1),linextrap);
	%test(ix)   = sqrt(trapz(uu,(uu-g(:,2*ix)).^2));
    end
  end
else % multip==0
  if (lower(def(1))=='n'||lower(def(1))=='e' || nargout>2),
    tp=[];mM=[];
    for ix=1:(m-1),
      tmp=dat2tp(xx(:,[1 ix+1]));
      tp=[tp; tmp];
      mM=[mM;tp2mm(tmp)];
    end
    
    cross1 = mm2lc(mM,cdef,0); %want upcrossings and maxima  
    cmax   = max(cross1(:,2));
    irr    = length(mM)/cmax;% approximately Tz/Tmaxima
  end
  if (lower(def(1))=='m')
    Fx  = edf(xx(n+1:end),'wdata',false);
    if plotflag==2
      stairs(Fx(:,1),Fx(:,2))
      pause(opt.delay)
    end
  end
  switch lower(def(1:3))
    case {'her'},
      ga1 = skew(xx(n+1:end));
      ga2 = kurt(xx(n+1:end))-3;
      up  = min(4*(4*ga1/3).^2,13);
      lo  = ga1^2*3/2;
      kurt1 = min(up,max(ga2,lo))+3;
      phat = [sa ga1,kurt1,ma ];
      [g, test] = hermitetr(g(:,1) ,phat);
      g2=[];
    case {'och'}
      phat = [sa skew(xx(n+1:end)) ma];
      [g test] = ochitr(g(:,1) ,phat);
      g2=[];
    case {'non'},
      [g, test, g2] = lc2tr(cross1,ma,sa,opt);%csm1,csm2,param,[plotflag monitor]);
    case {'emp'},  
      [g2, test, g] = lc2tr(cross1,ma,sa,opt);%csm1,csm2,param,[plotflag monitor]);
      test = sqrt(trapz(uu,(uu-smooth((g(:,1)-ma)/sa,g(:,2),1,uu,1)).^2));
    case {'mno'},
      [g, test, g2] = cdf2tr(Fx,ma,sa,opt);
    case  {'mem'},
      [g2, test, g] = cdf2tr(Fx,ma,sa,opt);
      test = sqrt(trapz(uu,(uu-smooth((g(:,1)-ma)/sa,g(:,2),1,uu,1)).^2));
  end 
end
