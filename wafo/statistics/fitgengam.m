function [phat]=fitgengam(data, varargin)
%FITGENGAM Parameter estimates for Generalized Gamma data.
%
% CALL:  [phat] = fitgengam(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ls'  = Least squares on log scale 
%                'ml'  = Maximum Likelihood method (default)
%                'mps' = Maximum Product of Spacings method. (robust)
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%                  
% 
% Example:
%   R = rndgengam(2,2,2,1,100);
%   phat = fitgengam(R);
%   plotfitsumry(phat);
%
%   close all;
%
% See also   pdfgengam,  cdfgengam, invgengam, rndgengam, momgengam

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


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on: Matlab 5.3
% History:
% revised pab Aug 2007
% -added fitgengamb as a subfunction
% revised pab 23.01.2004
%  - added chat0 
% revised pab 27.01.2001
%  - Improved the stability of estimation by adding a linear search to
%    find the largest sign change with positive derivative.
% revised pab 13.11.2000
%  - added check on estimated parameters -> return NaN when any is less
%    than zero
% revised pab 
%  -changed the order of a,b,c
%  -safer call to fzero 
% rewritten ms 10.08.2000

global WAFO_WSTATS_DEFAULT_PLOTFLAG
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'monitor',false, 'start',[],... 
  'optimset',optimset('disp','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',500)); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);


if any(data<=0)
  error('data must be strictly positive!')
end


data  = sort(data(:)); % make sure it is a vector

% Use Weibull start for c (assuming a=1)
logdata   = log(data);
meanlogdata = mean(logdata);

if ~isempty(options.start)
  start = options.start;
else
  % Do linear search to find a sign change pab 28.01.2001
  % -> more stable estimation
  
  start = 1./(6^(1/2)/pi*std(logdata)); % approx sqrt(a)*c
  
  c  = [linspace(max(eps,start/100),start,25) ...
    linspace(start+1,10*start+20,15)];
  
  [L,a] = fitgengamml(c,data,logdata,meanlogdata,options);
  
  ind = isnan(L)|(a<=0);
  if any(ind)
    L(ind) = [];
    c(ind)  = [];
  end
  dL=diff(L)./diff(c);
  %figure(1),  plot(c,L,(c(1:end-1)+c(2:end))/2 ,dL,'r');axis([0 inf -0.1 0.1]),figure(2)
  sl = sign(L);
  %choose the largest sign change with positive derivative
  ind = find(sl(1:end-1)<sl(2:end) , 1, 'last' );
  if ~isempty(ind)
    start = c(ind:ind+1);
  elseif any(dL>0)
    start = c(find(dL>0, 1 ));
  end
end

if strcmpi(options.method,'ml'),  % Maximum Likelihood
  chat = fzero(@(x)fitgengamml(x,data,logdata,meanlogdata,options),start,options.optimset);
  [L,ahat,bhat] = fitgengamml(chat,data,logdata,meanlogdata,options);
  phat0 = [ahat,bhat,chat];
elseif  strcmpi(options.method,'mps'),  % Maximum product spacing
  
  chat = mean(start);
  [L,ahat,bhat] = fitgengamml(chat,data,logdata,meanlogdata,options);
  phat0 = [ahat,bhat,chat]; 
  phat0 = fminsearch(@(x)logps(x,data,@cdfgengam),phat0,options.optimset);
elseif strcmpi(options.method,'ls')
   N = length(data);
   F = (0.5:N-0.5).'/N0;
   logR = log1p(-F);
   chat = start;
   [L,ahat,bhat] = fitgengamml(chat,data,logdata,meanlogdata,options);
   phat0 = [ahat,bhat,chat];
   phat0 = fminsearch(@(x)fitgengamls(x,data,logR,options),phat0,options.optimset);
else
  error(['Unknown method ' options.method '.']);
  %ahat = ahat0;  
end
  
if any(phat0<=0|isnan(phat0)), 
  phat0(1:3)=NaN;
  pcov=phat0;
  ciL=pcov;
  ciU = pcov;
else
  ahat = phat0(1);
  bhat = phat0(2);
  chat = phat0(3);
  
  c1 = psi(0,ahat);
  c2 = psi(0,ahat+1);
  c3 = psi(1,ahat);
  c4 = psi(1,ahat+1);
  pcov=[c3, -c1/chat, chat/bhat;
  -c1/chat, (1+ahat*c4+ahat*c2^2)/chat^2, -(1+ahat*c1)/bhat;
  chat/bhat, -(1+ahat*c1)/bhat, ahat*chat^2/bhat^2]/length(data);
   

 zcrit = -invnorm(options.alpha/2);
 pvar = diag(pcov).';
 ciL = phat0-sqrt(pvar).*zcrit;
 ciU = phat0+sqrt(pvar).*zcrit; 
  
  
end

[LPS,pvalue] = logps(phat0,data,'cdfgengam');
  phat = createfdata(options,'dist','pdfgengam',...
    'params',phat0,'lower',ciL,'upper',ciU,...
    'covar',pcov,'var',diag(pcov).',...
    'dataname',inputname(1),'data',data,...
    'loglikemax', -loglike(phat0,data,'pdfgengam'),...
    'logpsmax',-LPS,'pvalue',pvalue,'note',sprintf('Moran''s statistic on fit: pvalue = %g',pvalue));
  



if options.plotflag
  plotfitsumry(phat,options.plotflag)
end

phat = fdata(phat);



  function [L] = fitgengamls(phat,x,logR,options)
 % x   = data;
 tmp = sqrt(-logR);
 tmp2 = sqrt(-log1p(-cdfgengam(x,phat(1),phat(2),phat(3))));
 if options.monitor
   plot(x,[ tmp tmp2]); drawnow
 end,
 df = abs(tmp-tmp2).^(2);
 df(~isfinite(df)) = log(realmax);
 L = mean(df);
 
function [L,a1,b] = fitgengamml(c,data,logdata,meanlogdata,options)
% FITGENGAMML Maximum likelihood
    
    %ld  = log(data); 
    %mld = mean(logdata);
    sdc  = zeros(size(c));
    sdcld = sdc;
    
    for ix=1:numel(c)
      dc  = data.^c(ix);
      sdc(ix) = sum(dc);
      sdcld(ix) = sum(dc.*logdata);
    end
    a   = -1./(c.*(meanlogdata-sdcld./sdc));
    if nargout>1
      a1 = a;
    end
    a(a<=0) = nan; % Avoid error with wpsi/ gammaln for a<0 pab 27.01.2001
    n = numel(data);
    L = log(n*a) - psi(a) + c*meanlogdata -log(sdc);
    if nargout>2
      b = (sdc./n./a).^(1./c);
    end
    if options.monitor
      disp(['err = ' num2str(L,10)   ' phat = ' num2str(c,4) ])
    end

% function l=fitgengamb(b,data,F,def)
% %FITGENGAMB Is an internal routine for fitgengam
% %
% 
% % History
% % revised pab 21.01.2004
% % revised pab 27.01.2001
% if nargin<4||isempty(def),def=1;end
% 
% monitor = false;
% switch def
%   case 1,
%     % MLE
%     
%     %ld  = log(data); 
%     ld  = F;
%     mld = mean(ld);
%     sdb  = zeros(size(b));
%     sdbld = sdb;
%     
%     for ix=1:numel(b)
%       db  = data.^b(ix);
%       sdb(ix) = sum(db);
%       sdbld(ix) = sum(db.*ld);
%     end
%     a   = -1./(b.*(mld-sdbld/sdb));
%     a(a<=0) = nan; % Avoid error with gammaln for a<0 pab 27.01.2001
%     n = numel(data);
%     l = log(a) - psi(0,a) + b*mld -log(sdb/n);
%    
%   case 2, % LS-fit to empirical CDF
%     %LS
%     x   = data;
%     tmp = sqrt(-log1p(-F));
%     tmp2 = sqrt(-log1p(-cdfgengam(x,b(1),b(2),b(3))));
%     if monitor
%       plot(x,[ tmp tmp2]); drawnow
%     end
%     l = mean(abs(tmp-tmp2).^(2));
%   case 3,% Moment fit: data = E(x^2)/E(x)^2,F= E(x^3)/E(x^2)^(3/2)
%     % MOM
%     l = ...
% 	sum((gamma(b(1))*gamma(b(1)+2/b(2))/gamma(b(1)+1/b(2))^2-data)^2+...
% 	  (sqrt(gamma(b(1)))*gamma(b(1)+3/b(2))/gamma(b(1)+2/b(2))^1.5-F)^2);
%   case 4,% Moment fit: data = E(x^3)/E(x^2)^(3/2),F= E(x^4)/E(x^2)^(2)
%     l = sum((sqrt(gamma(b(1)))*gamma(b(1)+3/b(2))/gamma(b(1)+2/b(2))^1.5-data)^2+...
% 	  (gamma(b(1))*gamma(b(1)+4/b(2))/gamma(b(1)+2/b(2))^2-F)^2);
% end
% 
% if monitor
%   disp(['err = ' num2str(l,10)   ' phat = ' num2str(b,4) ])
% end




