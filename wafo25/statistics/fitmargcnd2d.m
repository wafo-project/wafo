function  [phat] =fitmargcnd2d(V,H,dist,res,method,monitor,chat0)
%FITMARGCND2D Parameter estimates for DIST2D data.
%
% CALL:  phat = fitmargcnd2d(x1,x2,dist,[res,noverlap C],method,monitor);
%
%        phat = structure array containing
%               x      = cellarray of distribution parameters
%               CI     = Confidence interval of distribution parameters
%               dist   = as explained below
%               method = ------||----------
%               note   = string
%               date   = date and time of creation 
%       x1,x2 = input data
%        dist = list of distributions used in the fitting of X1 given X2 
%               and X2, respectively. Options are 'tgumbel', 'gumbel', 
%               'lognormal','rayleigh','trayleigh','weibull','tweibull',
%               'gamma','ggamma'.
%         res = resolution in the conditional fitting  (default range(x2)/12)
%    noverlap = 0,1,2,3..., nr. of groups above and below to overlap in
%               the conditional fitting (default 0)
%           C = defines the location parameter by the following:
%               Cv = min({X1|X2=x}-C,0)
%               if no value is given then Cv=0. 
%      method = list containing method used for fitting of 
%               X1 given X2 and X2, respectively.
%               Options 'MLE' or 'MOM'  (default {'MLE',MLE'})
%     monitor = 1 if monitor the fitting  (default=0) 
%
%  FITMARGCND2D fits the following distribution: F{x1|X2=x2}*F{x2}. 
%  The parameter(s) of the unconditional distribution of X2,
%  F{x2}, is  returned in phat.x{2}. The parameters of the conditional
%  distribution of X1 given X2 are returned in phat.x{1}. The first column
%  in PHAT.X{1} contains the X2 values the parameters in column 2 and 3 are
%  conditioned on.  PHAT.CI{1} and PHAT.CI{2} gives 95% confidence
%  interval if MLE are selected for METHOD{1} and METHOD{2}, respectively.
%
% Example:
%   R = rndray(2,1000,2); x1 = linspace(0,10)';
%   phat0=struct('x',[]); phat0.x={[x1,2*ones(size(x1))] 2 };
%   phat0.dist={'rayl','rayl'};
%   phat = fitmargcnd2d(R(:,1),R(:,2),{'rayl','rayl'}); 
%   sphat = margcnd2dsmfun2(phat,x1,0); % smooth the parameters
%   plotmargcnd2dfit(phat,sphat); 
%   figure(2) % compare fitted model with theoretical
%   f = pdfmargcnd2d2(x1,x1,phat); 
%   fs = pdfmargcnd2d2(x1,x1,sphat); 
%   fe = pdfmargcnd2d2(x1,x1,phat0); 
%   pdfplot(fs); hold on,
%   pdfplot(fe,'k--')
%   plot(R(:,1),R(:,2),'.'),hold off
%
% See also rndmargcnd2d, pdfmargcnd2d, cdfmargcnd2d prbmargcnd2d

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% tested on: matlab 5.2
% history:
% revised pab Sept 2005
% -fixed abug: [a{1}] = 1 is not allowed anymore!! replaced with a{1} = 1;
% revised pab Jul2004
% Added secret option of C0 to fitgengam 
% revised pab 03.12.2000
% -added truncated weibull and rayleigh
% revised pab 12.11.2000
%  - added pdfgengam option
%  - For monitor> 0:  removed normplot ...etc and replaced it with a call
%     to edf instead
% Per A. Brodtkorb 20.10.1998

error(nargchk(3,7,nargin))
ptime     = 1; %pause length if monitor=1
printflag = 0;%print if monitor=1

[HI I] = sort(H(:));%sorting H
VI     = V(I);

if (nargin< 3)||isempty(dist),
  dist={'tgumbel','rayleigh'};
else
  nd=length(dist);
  if nd==1,
  dist={dist{1},'rayleigh'};
  end
end

if (nargin<5)||isempty(method)
  method={'MLE','MLE'}; % options MLE or MOM
else
  nd=length(method);
  if nd==1,
    method={method{1}, 'MLE'};
  end
end

if (nargin<6)||isempty(monitor)
  monitor=0; %if 1 monitor the fitting 
end
%monitor=1;


noverlap=0;C=[];
if nargin<4||isempty(res)
  res=range(H(:))/12; % resolution
else
  if length(res)>1
    noverlap=res(2);
  end
  if length(res)>2
    C=res(3);
  end
  res=res(1);
end

%N=length(HI);

grp=floor(HI/res)+1; % dividing the data into groups

if monitor
  [len] = bincount(grp);
  utp = 1/(2*max(len));
  ax = [0 max(VI),utp,1-utp];
  name1=inputname(1);
  name2=inputname(2);
  if isempty(name1), name1 = 'x1'; end
  if isempty(name2), name2 = 'x2'; end
  
end
[ phhat, phci]=distfit(HI,dist{2},method{2},monitor);%unconditional fitting
if monitor
  xlabel(name2)
  pause(ptime)
  if printflag, 
    %print('-Pmhlaser'); 
  end  
end

Ngrp=grp(end); % # groups

if strcmpi(dist{1}(1:2),'ra'),
  npar=1;
elseif strcmpi(dist{1}(1:2),'gg')||strcmpi(dist{1}(1:2),'tw'),
  npar = 3;
else
  npar=2;  
end
%C = input('Specify location parameter (default none)'); 
if isempty(C),  
  pvhat=zeros(Ngrp,1+npar);
else
  Cvold=0;
  pvhat=zeros(Ngrp,2+npar);
end

pvci=zeros(Ngrp,2*npar);

pvhat(:,1)=linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
if nargin<7
  chat0 = [];
end

Nmin = min(max(6,Ngrp),25); % Minimum number of data in groups
m = zeros(Ngrp,1);
v = m;
Ni = m;
for ix=1:Ngrp,
  ind  = (grp==ix);  
  tmp  = VI(ind);
  Ni(ix) = length(tmp);
  if length(tmp)>max(4,0),% if less than 4 observations in the group 
    m(ix)=mean(tmp); % mean of data in group ix
    v(ix)=std(tmp).^2;  % variance of data in group ix
  else 
    m(ix)=NaN;
    v(ix)=NaN;
  end
  %if 1 | ((ix-1)*res>2)
  %  chat0 = 1.5;
  %end
  switch class(chat0)
   case 'inline'
    chat00 = chat0(pvhat(ix,1));
   case 'struct'
    chat00 = ppval(chat0,pvhat(ix,1));
   otherwise
    chat00 = chat0; 
  end
  
  if ix>=Ngrp-noverlap+1 || ix<=noverlap
    for iy=1:min(min(noverlap,ix-1),Ngrp-ix)
      kup = (grp==ix+iy);
      if any(kup)
	ind = (ind | (grp==ix-iy) | kup);
      end
    end 
  end             
  
  tmp=VI(ind);%find data in group number ix
  
  if length(tmp)<Nmin,% if less than Nmin observations in the group 
    %grp(grp==ix)=ix-1; %merge the data with group below
    %if (ix>3)&~isempty(tmp), grp(grp==ix-2)=ix-1;end % also merge the
                                                     % data in grp ix-2 with grp ix-1
    pvhat(ix,:)=NaN;
    pvci(ix,:)=NaN;
  else
    
    if ~isempty(C)
      %disp('Cv')
      Cv=max(min(tmp)-C,0);
      
      Cv=min(Cv,Cvold+C/30);
      tmp=tmp-Cv;
      Cvold=Cv;
      pvhat(ix,end)=Cv;
    end
    % conditional fitting 
    [ pvhat(ix, 2:npar+1),...
      pvci(ix,1:2*npar)]=distfit(tmp(:),dist{1},method{1},monitor,chat00);
    if monitor
      axis(ax)
      xlabel([name1 ' |  '  name2 '=' num2str(pvhat(ix,1)) ])
      pause(ptime)
      if printflag, 
        %print -Pmhlaser ; 
      end   %print -Pmhlaser
    end
  end
end

%phat=struct('x',[],'CI',[],'dist',dist,'method',method);
phat.x=cell(2,1);
phat.CI=cell(2,1);
phat.x{1}=pvhat;
phat.x{2}=phhat;
phat.CI{1}=pvci;
phat.CI{2}=phci;
phat.dist=dist;
phat.method=method;
phat.note=[];
phat.date=datestr(now);
phat.res=res;
phat.noverlap=noverlap;
phat.C=C;
phat.visual=[];
phat.csm=[];
phat.lin=[];
phat.stats1{1} = pvhat(:,1); % x2
phat.stats1{2} = m; %conditional mean given x2
phat.stats1{3} = v; % conditional variance given x2
phat.stats1{4} = Ni;% conditional number of samples given x2
%pvhat(end,:)=[];
%pvhat(1,:)=[];
%pvci(end,:)=[];
%pvci(1,:)=[];
function [pvhat,pvci]=distfit(tmp,dist,method,monitor,chat)
  if nargin>4 && ~isempty(chat) && isfinite(chat)
    chat0 = chat;
  else
    chat0 = [];
  end
  if strncmpi(method,'mle',3),
    switch lower(dist(1:2)),
      case 'tr', [pvhat, tmp2] = fitraymod(tmp,monitor) ;
      case 'ra', [pvhat, tmp2] = fitray(tmp,monitor) ;
      case 'tg', [pvhat, tmp2] = fitgumb(tmp,monitor);
      case 'gu', [pvhat, tmp2] = fitgumb(tmp,monitor);
      case 'lo', [pvhat, tmp2] = fitnorm(log(tmp),monitor);
      case 'ga', [pvhat, tmp2] = fitgam(tmp,monitor);
     case 'gg',  [pvhat, tmp2] = fitgengam(tmp,monitor,chat0);
	%if any(pvhat<0),pvhat(1:3) =NaN; end
      case 'we', [pvhat, tmp2] = fitweib(tmp,monitor);
      case 'tw', [pvhat, tmp2] = fitweibmod(tmp,monitor);
	%[tmp tmp2]=gumbfit(log(tmp));
	%pvhat=[exp(tmp(2)) 1/tmp(1) ]
      otherwise, error('Unknown distribution')
    end
    if any(isnan(pvhat))
      pvci=[pvhat pvhat];
    else
      tz=size(tmp2);
      alpha=0.05; % 95% confidense interval
      pint = [alpha/2; 1-alpha/2];
      if tz(1)==tz(2),
        sa = diag(tmp2)';
      else
        sa  = tmp2(:)';
      end
      nt=length(sa);
      %pvhat
      pvci = invnorm(repmat(pint,[1 nt]),[pvhat;pvhat],[sa;sa]);
      pvci=pvci(:).';
    end
    if 0 %monitor
      switch lower(dist(1:2)),
	case 'ra',plotray(tmp);
	case 'gu',plotgumb(tmp);
	case 'tg',plotgumb(tmp);
	case 'lo',plotnorm(log(tmp));
	case 'ga',plotresq(tmp,'pdfgam');
	case 'gg',plotresq(tmp,'pdfgengam');
	case 'we',plotweib(tmp);
	otherwise, error('Unknown distribution')
      end    
    end
    
  else  % MOM fit
   
    switch lower(dist(1:2))
     case {'tg','gu'} ,  pvhat =plotgumb(tmp);
      
     case 'we', pvhat = plotweib(tmp);
      
     case 'lo', pvhat = plotnorm(log(tmp));
      
     case 'ga', ma=mean(tmp);sa=var(tmp);
       pvhat=[ma^2/sa sa/ma];
      
     case 'ra', error('MOM is not implemented for Rayleigh distribution')	 
     otherwise , error('Unknown distribution')
    end
    pvci = repmat(nan,1,2*length(pvhat));
    
    if monitor
      switch lower(dist(1:2)),
	case {'gu' 'tg'}, plotgumb(tmp);
	case 'lo',  plotnorm(log(tmp));
	case 'ga', plotresq(tmp,'pdfgam');
	case 'we',  plotweib(tmp);
	otherwise, error('Unknown distribution')
      end
    end
  end
  return





