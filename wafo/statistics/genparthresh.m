function [u,shape,scale,pvalue] = genparthresh(data,varargin)
%GENPARTHRESH Automatic Threshold selection for GPD data
%
% CALL  [u,shape,scale,pvalue] = genparthresh(data,options)
%
% u       = selected threshold
% shape, 
%   scale = WDATA objects with estimated GPD parameters as function of
%           threshold.
% data    = vector of data.
% options = options structure defining the range of thresholds.
%        .method : Method used in the fitting
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 10).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),20))
%        .alpha: Confidence coefficient (default 0.05)
%
% GENPARTHRESH estimate GPD model parameters over a range of thresholds and
% automatically selects the most appropriate one.
% The purpose is to determine the threshold where the upper tail of the data 
% can be approximated with the generalized Pareto distribution (GPD). The 
% GPD is appropriate for the tail, if the shape- and modified scale- 
% parameter is constant.
%
% Example
%  opt = genparthresh('defaults');
%  opt.Nu = 20;
%  R = rndgenpar(0.1,2,2,100,1);
%  [u,shape,scale,pvalue] = genparthresh(R-2,opt); 
%  plot(shape), figure(gcf+1)
%  plot(scale), figure(gcf+1)
%  plot(pvalue)
%
%  close all;
%
% See also fitgenpar, reslife, disprsnidx

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


%disp('genparthresh')
u = [];
options = struct('method','MPS','shapefix',nan,'shapemin',0.001,'umin',[],'umax',[],'Nu',10,'Nmin',10,'alpha',0.05,'plotflag',0);

if ischar(data) && strcmpi(data,'defaults')
  shape = options;
  return
end
options = parseoptions(options,varargin{:});

sd = sort(data);
n = length(data);

if isempty(options.Nmin)
   options.Nmin=10;
else
  options.Nmin = max(options.Nmin,0);
end
if options.Nmin>n/2
  warning('WAFO:FITGENPARRANGE','Nmin too large!')
  options.Nmin = min(options.Nmin,ceil(n/2));
end
try
  sdmax = sd(max(n-options.Nmin,5));
catch
  sdmax = sd(max(n-1,1));
end
if isempty(options.umax)
  options.umax = sdmax;
else
  options.umax = min(options.umax,sdmax);
end
xMedian = sd(max(floor(n/2),1));
if isempty(options.umin)
  options.umin = xMedian;
else
  options.umin = max(options.umin,sd(1));
end


if isempty(options.Nu)
  options.Nu = min(n-options.Nmin,20);
end
shapefix = options.shapefix;
if ~isnan(shapefix)
  options.shapemin = shapefix;
end

u1 = linspace(options.umin,options.umax,options.Nu).';


nan1 = nan;
pvalue1 = nan1(ones(options.Nu,1));
phat1 = nan1(ones(options.Nu,2));
s = phat1;
num = ones(options.Nu,1);
method = options.method;

%xmax = sd(n);
%returnPrbs = fliplr([1/(n*200) logspace(-log10(n*20),log10(.5))]);
%RL = invgenpar(returnPrbs/n,0.1,2,2,'lowertail',false);
%semilogx(1./returnPrbs,RL,'r.'), hold on

returnPrbs = [1/(n*20),1/(n*200)];

linscale = true;
max_point = zeros(1,options.Nu);
tmp = zeros(options.Nu,length(returnPrbs));
for ix=1:options.Nu;
  [phat] = fitgenpar(data(data>u1(ix)),'method',method,'fixpar',[shapefix,nan,u1(ix)]);
   phat1(ix,:) = phat.params(1:2);
   pvalue1(ix) = phat.pvalue;
   pcov = phat.covariance(1:2,1:2);
   phat1(ix, 2) = phat1(ix, 2) + phat1(ix, 1) * u1(ix); % modified scale
   if linscale % linear scale
     d = [u1(ix);1];
   else
     % Assume gaussian on logarithmic scale.
     d = [u1(ix);1]/phat1(ix,2);
     phat1(ix,2) = log(phat1(ix,2));
   end
   s(ix,:) = sqrt(diag(pcov).');
   s(ix,2) = sqrt(d.' * pcov *d);
   num(ix) = sum(data>u1(ix));
  if options.shapemin<=phat.params(1)
     
    [returnLevels] = invgenpar(returnPrbs/num(ix),phat,'lowertail',false);
    q1000 = returnLevels(end-1);
    q10000 = returnLevels(end);
    qdiff = q10000-q1000;
    
    if ~isnan(options.shapefix) || abs(qdiff) < max(0.2*sqrt(50/n),0.12)*abs(q10000)
        max_point(ix) = phat.params(2)/phat.params(1)+phat.params(3);
        tmp(ix,:) = returnLevels;
        
        %semilogx(1./returnPrbs,returnLevels,'r');%,1./returnPrbs,qlo,'b',1./returnPrbs,qup,'b'), hold on
        %hold('on')
       % shg 
        %u(ix)
        %pause(1)
        
    end
  end
end
%max_point

if any(tmp)
  dim = 1;
  if isnan(options.shapefix)
    mRL = sum(tmp,dim)./sum(tmp~=0,dim);
  else
    mRL = max(tmp,[],dim);
  end
  [df,ix] =  min(abs(tmp(:,end-1)-1.15*mRL(end-1)));
  [df1,ix1] =  min(abs(tmp(:,end)-1.15*mRL(end)));
    %semilogx(1./returnPrbs,mRL,'ro'), hold on
    
    u = max(u1(ix),u1(ix1));
end
if nargout>1 || options.plotflag>0
p = 1-options.alpha;
alpha2 = options.alpha/2;
% Approximate P% confidence interval
if 0,
 Za = invnorm(alpha2);   % assume known variance
 
 up = phat1 + Za.*s;
 lo = phat1 - Za.*s;


else
  Za   = invt(alpha2,num-1); % unknown variance
  up = phat1 + Za(:,ones(1,2)).*s;
  lo = phat1 - Za(:,ones(1,2)).*s;
  
  
end
if ~linscale
  phat1(:,2) = exp(phat1(:,2));
  up(:,2)  = exp(up(:,2));
  lo(:,2)  = exp(lo(:,2));
end



titleTxt1 = sprintf('GPD shape parameter with %d%s CI',100*p,'%');
titleTxt2 = sprintf('GPD modified scale parameter with %d%s CI',100*p,'%');
titleTxt3 = sprintf('P-value of fitted GPD');

shape = createwdata('data',phat1(:,1),'args',u1,...
'dataCI',[lo(:,1),up(:,1)],'title',titleTxt1,'labels',{'Threshold','Shape'},...
  'workspace',options,'note',titleTxt1);


scale = createwdata('data',phat1(:,2),'args',u1,...
'dataCI',[lo(:,2),up(:,2)],'title',titleTxt2,'labels',{'Threshold','Modified Scale'},...
  'workspace',options,'note',titleTxt2);

pvalue = createwdata('data',pvalue1,'args',u1,...
'title',titleTxt3,'labels',{'Threshold','P-value'},...
  'workspace',options,'note',titleTxt3);


  shape = wdata(shape);
  scale = wdata(scale);
  pvalue = wdata(pvalue);
  if options.plotflag>0
    subplot(2,1,1)
    plot(shape,options.plotflag,'.')
    subplot(2,1,2)
    plot(scale,options.plotflag,'.')
  end

end


