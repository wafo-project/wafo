function [phat]=fitweib(data,varargin)
%FITWEIB Parameter estimates for Weibull data.
%
% CALL:  phat = fitweib(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ml'  = Maximum Likelihood method (default)
%                'mps' = Maximum Product Spacing method
%                'ls'  = Least squares on log(1-F) scale
%     .fixpar  : vector giving the fixed parameters. (Must be empty or 
%               have the same length as the number of parameters. 
%               Non-fixed parameters must then be given as NaN's)
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
%   R=rndweib(10,2,1,sz);
%   [phat] = fitweib(R)
%   x = linspace(0,25,200).';
%   [p,plo,pup] = cdfweib(x,phat);
%   plotflag = 1012
%   plotedf(R,[],plotflag), hold on
%   f = wdata(p,x); set(f,'dataCI',[plo pup])
%   plot(f,'r',plotflag-1), hold off
%
% See also  cdfweib, pdfweib, invweib, rndweib, momweib

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
% and Life Span Models", p. 25 ff, Marcel Dekker.

%Tested on: matlab  5.3
% History:
% revised pab 03.11.2000
% - added
% revised pab 24.10.2000
%  - added nargchk + safer call to fzero
%  - made sure data is a vector 
% rewritten ms 20.06.2000

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','fixpar',[nan, nan, 0],'alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'monitor',false,'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

data  = data(:);                            % make sure it is a vector

somefixed = ~isempty(options.fixpar);
if somefixed
  phat2 = options.fixpar;
  isnotfixed = (~isfinite(phat2));
  i_fixed   = find(isfinite(phat2));
  if isfinite(phat2(3))
    c0 = phat2(3);
  else
    c0 =  min(data)-0.01*std(data);
  end
else
  c0 = min(data)-0.01*std(data);
end
[phat1] = fitml2(data-c0,options);
phat0 =  [phat1,c0];
if somefixed
  phat0 = phat0(isnotfixed);
end
if strcmpi(options.method,'ml'),  % Maximum Likelihood
  dosearch = somefixed && any(i_fixed<3);
elseif  strcmpi(options.method,'mps'),  % Maximum product spacing
  dosearch = true;
elseif strcmpi(options.method,'ls'),
  dosearch = false;
  phat0 = fitls(data,phat0,options);
else
  error(['Unknown method ' options.method '.']);
end

phat = mlest(@pdfweib,phat0,data,{options,'search',dosearch});
phat.dataname = inputname(1);


function [phat,pcov] = fitml2(data,options)
start = 1./(6^(1/2)/pi*std(log(data)));


chat = fzero(@fitcweib,start,options.optimset,data);


ahat = mean(data.^chat).^(1./chat);
phat = [ahat(:), chat(:)];
if nargout>1
  pcov=[1.109*ahat^2/chat^2,0.257*ahat;0.257*ahat,0.608*chat^2]/length(data);
end

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

def=3; %  What is def? See 'weibfun' 
fixedPhat = options.fixpar;
if isempty(fixedPhat)
 fixedPhat = [nan nan nan];
end
phat = fminsearch(@weibfun,phat0,options.optimset,sd,F,def,monitor,fixedPhat);





function y=weibfun(phat,x,F,def,monitor,fixedPhat)
% WEIBFUN Is an internal routine for fitls
%

fixedPhat(~isfinitefixedPhat) = phat;

a =fixedPhat(1);
b =fixedPhat(2);
c1 =fixedPhat(3);

c = -(c1);
N = length(F);
xmin = min(x(:));

penalty = abs(N*c*(xmin<c1));

%monitor = logical(1);
switch def
  case 1, % fit to sqrt(-log(1-F))
  if monitor
    plot(x,sqrt(-log(1-F)),....
	x,sqrt(((x+c)./a).^b)); drawnow
  end
  
  y=mean((-sqrt(-log(1-F))+...
      sqrt(((x+c)./a).^b).^2))*(1+penalty) + penalty;
case 2, % fit to (-log(1-F))
  if monitor
    plot(x,(-log(1-F)),...
	x,(((x+c)./a).^b)); drawnow
  end
  
  y=mean((-(-log(1-F))+...
      (abs((x+c)./a).^b)).^2)*(1+penalty) + penalty;

case   3, % fit to (-log(1-F)).^(1/b)
  if monitor
    plot(x,(-log1p(-F)).^(1/b),x,...
	(((x+c)./a))); drawnow
  end
  y=mean((-(-log(1-F)).^(1/b)+...
      (((x+c)./a))).^2)*(1+penalty)+penalty;
case 4,  % fit to (-log(1-F)+(c/a).^b).^(1/b)
  if monitor
    plot(x,(-log(1-F)).^(1/b),x,(x+c)./a); drawnow
  end
  y=mean((-(-log(1-F)).^(1/b)+(x+c)./a).^2)*(1+penalty)+penalty;

end

if monitor
  disp(['err = ' num2str(y,10)   ' a b c = ' num2str([a,b,c1],4) ])
  shg
end

function l=fitcweib(c,data)
%FITCWEIB Is an internal routine for fitweib
%

l=1./c-sum(data.^c.*log(data))/(sum(data.^c))+sum(log(data))/size(data,1);
