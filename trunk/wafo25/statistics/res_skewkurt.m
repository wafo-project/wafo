function mrl=res_skewkurt(data,varargin)
%RES_SKEWKURT Residual L-kurtosis to L-skewness ratio  vs thresholds
%
% CALL  lmr = res_skewkurt(data,options)
%
% lmr    = L-moment ratio for excesses over thresholds, u.
%
% data    = vector of data.
% options = options structure defining estimation of mrl. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 3).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),100))
%        .alpha: Confidence coefficient (default 0.05)
%        .plotflag:
%
% RES_SKEWKURT estimate L-moment ratio for excesses over thresholds. The 
% purpose of LMR is to determine the threshold where the upper tail of the 
% data can be approximated with the generalized Pareto distribution (GPD).
% The GPD is appropriate for the tail, if the LMR is one as function of
% threshold, u. Theoretically in the GPD model the ratio
%
%   KU*(5+SK)/(SK*(1+5*SK)) = 1
%
% where KU and SK are L-kurtosis and L-skewness, respectively.
%
% Example
%  opt = res_skewkurt('defaults'); % Get default options
%  opt.Nu = 20;
%  xn = load('sea.dat');
%  Ie = findpot(xn,0,5);
%  mrl = res_skewkurt(xn(Ie,2),opt);
%  plot(mrl)
%  
% See also reslife, fitgenpar, fitgenparrange, disprsnidx


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



options = struct('umin',[],'umax',[],'Nu',100,'Nmin',5,'Nboot',2000,'alpha',0.05,'plotflag',0);

if ischar(data) && strcmpi(data,'defaults')
  mrl = options;
  return
end
options = parseoptions(options,varargin{:});

sd = sort(data);
n = length(data);

if isempty(options.Nmin)
   options.Nmin=5;
else
  options.Nmin = max(options.Nmin,0);
end
if options.Nmin>n/2
  warning('WAFO:RES_SKEWKURT','Nmin possibly too large!')
end

sdmax = sd(n-options.Nmin);
if isempty(options.umax)
  options.umax = sdmax;
else
  options.umax = min(options.umax,sdmax);
end

if isempty(options.umin)
  options.umin = sd(1);
else
  options.umin = max(options.umin,sd(1));
end


if isempty(options.Nu)
  options.Nu = min(n-options.Nmin,100);
end


u = linspace(options.umin,options.umax,options.Nu).';


nan1 = nan;
lmr = nan1(ones(options.Nu,1));
num = lmr;
lmrCI = lmr(:,ones(1,2));


for ix=1:options.Nu;
  [lmr(ix),lmrCI(ix,:),num(ix)] = myLMR(data(data>u(ix))-u(ix),options);
end
% 
options.CI = lmrCI;
options.numdata = num;
%
titleTxt = sprintf('Residual L-kurtosis to L-skewness ratio with %d%s CI',100*(1-options.alpha),'%');
mrl = createwdata('data',lmr,'args',u,...
    'dataCI',options.CI,...
'title',titleTxt,'labels',{'Threshold','Ratio'},...
    'workspace',options,'note',titleTxt);
  
if ~isoctave
  mrl = wdata(mrl);
  if options.plotflag>0
    plot(mrl,options.plotflag,'.')
  end
end



function [lmr,lmrCI,n] = myLMR(dat1,options)
% Return L-moment ratio ,i.e., ku*(5+sk)/(sk.*(1+5*sk))
% where ku = L-kurtosis and sk = L-skewness

x = dat1(:);
n = length(x);
[sk,ku] = lskewkurt(x);
lmr = ku*(5+sk)/(sk.*(1+5*sk));

% LMR confidence interval using bootstrap:
B = options.Nboot;
alpha  = options.alpha;

xB = zeros(n,B);
J = ceil(rand(n*B,1)*n);
xB(:) = x(J); 
   
[skB,kuB] = lskewkurt(xB);
lmrB = kuB.*(5+skB)./(skB.*(1+5*skB));

lmrCI = percentile(lmrB,[alpha/2,1-alpha/2]);
   
function [sk,ku] = lskewkurt(Q)

m = 4;
x = sort(Q);
n = size(x,1);

lmom = cell(1,m);
lmom{1} = mean(x);
  if m>1
    temp = ((1-n):2:(n-1));
    p = temp./n;
    lmom{2} = p*x/n;
    if m>2
       p1   = ones(1,n);
       tmp = p;
    end
  end

  for i = 3:m
    p2 = p1;
    p1 = p;
    p = ( (2*i-3)*tmp.*p1 - (i-2) * p2 ) /((i-1));
    %p = ( (2*i-3)*temp.*p1 - (i-2) * (n + i - 2)  * p2 ) /((i-1) * (n - i+1));
    lmom{i} = (p*x/n)./lmom{2}; 
  end
  sk = lmom{3};
  ku = lmom{4};
