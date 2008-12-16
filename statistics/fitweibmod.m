function phat =fitweibmod(data, varargin)
%FITWEIBMOD Parameter estimates for truncated Weibull data.
%
% CALL:  phat = fitweibmod(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ml'  = Maximum Likelihood method 
%                'mps' = Maximum Product Spacing method (default)
%                'ls'  = Least squares on log(1-F) scale 
%     .fixpar  : vector giving the fixed parameters. (Must be empty or 
%               have the same length as the number of parameters. 
%               Non fixed parameters must then be given as NaN's)
%               (default [nan nan 0])
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details))
% 
% Example:
%  sz = [1 100]
%   R=rndweib(2,2,0,sz);
%   Rt=R(R>1)-1;  % Truncated weibul with a=2, b=2, c=1
%   [phat] = fitweibmod(Rt)
%   plotfitsumry(phat)
%
% See also  cdfweibmod, pdfweibmod, invweibmod, rndweibmod, momweibmod

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 25 ff, Marcel Dekker.

%Tested on: matlab  5.3
% History:
% revised pab July2004  
% revised pab Dec2003
% revised pab 03.11.2000
% - added
% revised pab 24.10.2000
%  - added nargchk + safer call to fzero
%  - made sure data is a vector 
% rewritten ms 20.06.2000

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','fixpar',[nan, nan, nan],'alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'monitor',false,'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

data  = data(:);                            % make sure it is a vector


phat1 = fitweib(data);
phat0 = [phat1.params(1:2),0];
somefixed = ~isempty(options.fixpar);

if somefixed
  isnotfixed = (~isfinite(options.fixpar));
  phat0 = phat0(isnotfixed);
end
if strcmpi(options.method,'ml') || strcmpi(options.method,'mps'),
  % Maximum Likelihood or Maximum product spacing
  dosearch = true;
elseif strcmpi(options.method,'ls'),
  dosearch = false;
  phat0 = fitls(data,phat0,options);
else
  error(['Unknown method ' options.method '.']);
end

phat = mlest(@pdfweibmod,phat0,data,{options,'search',dosearch});
phat.dataname = inputname(1);



function [phat] = fitls(data,phat0,options)

N     = length(data);

sd = sort(data);

F = (0.5:1:(N - 0.5))'/N;
if N>10000
  Fi = linspace(F(1),F(end-5),10000).';
  %Fi = fliplr(logspace(log10(F(end)),log10(F(1)),10000)).';
  sd =interp1(F,sd,Fi,'linear');
  F  = Fi;
end

monitor = options.monitor;

def=2; %  What is def? See 'weibfun' 
fixedPhat = options.fixpar;
if isempty(fixedPhat)
 fixedPhat = [nan nan nan];
end
phat = fminsearch(@weibfun,phat0,options.optimset,sd,F,def,monitor,fixedPhat);





function y=weibfun(phat,x,F,def,monitor,fixedPhat)
% WEIBFUN Is an internal routine for fitls
%

fixedPhat(~isfinite(fixedPhat)) = phat;

a =fixedPhat(1);
b =fixedPhat(2);
c1 =fixedPhat(3);

c = abs(c1);

N = length(F);
%monitor = logical(1);
switch def
  case 1, % fit to sqrt(-log(1-F))
  if monitor
    plot(x,sqrt(-log(1-F)),....
	x,sqrt(((x+c)./a).^b-(c/a).^b)); drawnow
  end
  
  y=mean((-sqrt(-log(1-F))+...
      sqrt(((x+c)./a).^b-(c/a).^b)).^2) + N*c*(c1<0);
case 2, % fit to (-log(1-F))
  if monitor
    plot(x,(-log1p(-F)),...
	x,(((x+c)./a).^b-(c/a).^b)); drawnow
  end
  y=mean((-(-log(1-F))+...
      (((x+c)./a).^b-(c/a).^b)).^2)+N*c*(c1<0);

case   3, % fit to (-log(1-F)).^(1/b)
  if monitor
    plot(x,(-log(1-F)).^(1/b),x,...
	(((x+c)./a).^b-(c/a).^b).^(1/b)); drawnow
  end
  y=mean((-(-log(1-F)).^(1/b)+...
      (((x+c)./a).^b-(c/a).^b).^(1/b)).^2)+N*c*(c1<0);
case 4,  % fit to (-log(1-F)+(c/a).^b).^(1/b)
  if monitor
    plot(x,(-log(1-F)+(c/a).^b).^(1/b),x,(x+c)./a); drawnow
  end
  y=mean((-(-log(1-F)+(c/a).^b).^(1/b)+(x+c)./a).^2)+N*c*(c1<0);
case 5, % fit x/a to ((-log(1-F)+abs(a)).^(1/b));
       
  tmp = ((-log(1-F)+abs(a)).^(1/b))-abs(a)^(1/b);  
  p = (x\tmp).'; % Linear LS fit to find 1/a
  tmp = tmp/p(1);
  if monitor
    plot(x,x,x,tmp); drawnow
  end
  % Equal weigth on all x: 
  y = (mean(abs((tmp-x)).^(2)))+N*abs(a)*(a<0)+ (b-15)^2*(b>15)/N;
case 6, % fit x/a to ((-log(1-F)+abs(a)).^(1/b));
  
  cda = abs(a).^(1/b); % = c/a
  tmp = ((-log1p(-F)+abs(a)).^(1/b))-cda;
  p = (x\tmp).'; % Linear LS fit to find 1/a
  tmp = tmp/p(1);
  
  if 0 %monitor
    plot(x,x,x,tmp); drawnow
  end
  
  tmp3 =  (-log1p(-F));
  tmp4 = (((x*p(1)+cda)).^b-abs(a));
  
  % fit to (-log(1-F))  
  % More weight on the tails: The tail is fitted very well
  y = mean(abs(x-tmp).^(2)+abs(tmp3-tmp4).^(2))+N*abs(a)*(a<0)+(b-6)*(b>10)/N;
  if monitor
    plot(x,[x, tmp],x,[tmp3,tmp4]); drawnow
  end
case 7
  pac=[0.00077598974699  -0.02620368505187   1.28552709525102  -0.73037371897582];
  %pba=[-0.00641052386506   0.13245900299368   0.45810897559486  -0.38495820627853];
  % c = abs((a^1.25-0.4)/1.41)+.2;
  %c = abs((a^1.25-0.2)/1.45);
  %a = polyval(pba,b);
  c = polyval(pac,a);
  %c = abs((a^1.25-0.57)/1.41);
  cda = abs(c/a);
  tmp = (((-log(1-F)+cda^b).^(1/b))-cda)*a;
  %tmp3 =  (-log(1-F));
  %tmp4 = ((x+c)/a).^b-cda^b;
  if monitor
    plot(x,[x, tmp]); drawnow
  end
  y = mean(abs(x-tmp).^(2))+N*abs(a)*(a<=0)+(b-6)*(b>6)/N;
case 8
  tmp = sqrt(-log(1-F));
  tmp2 = sqrt(-log1p(-cdfweibmod(x,a,b,c)));
  if monitor
    plot(x,[ tmp tmp2]); drawnow
  end
  y = mean(abs(tmp-tmp2).^(2));
end

if monitor
  disp(['err = ' num2str(y,10)   ' a b c = ' num2str([a,b,c],4) ])
end

