function  [ei, Tc] = extremalidxrange(data,varargin)
%EXTREMALIDXRANGE Extremal Index vs threshold
%
%  CALL [EI,Tc] = extremalidxrange(data,options)
%
%  EI   = Extremal index
%  Tc   = minimum distance between clusters.
%  data = twocolumn matrix with sampled times in the first column
%         and values the second columns. Sampling frequency must be 1Hz.
% options = options structure defining estimation of EI. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 5).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),100))
%        .u    : (default linspace(umin,umax,Nu)
%        .alpha: Confidence coefficient (default 0.05)
%        .plotflag:
%
%   EXTREMALIDXRANGE estimate the Extremal Index (EI) which is a measure of
%   independence. The EI is one if the data are independent and less than
%   one if there are some dependence. The extremal index can also be
%   intepreted as the reciprocal of the mean cluster size.
%
%
% Example
%  xn = load('sea.dat');
%  umin = max(xn(:,2))/2
%  opt = extremalidxrange('defaults');
%  opt.umin = umin;
%  opt.umax = 1.6;
%  opt.Nu = 20;
%  Ie = findpot(xn,umin);
%  di = disprsnidx(xn(Ie,:),opt,'Tb', 100);
%  plot(di) % a threshold around 1 seems appropriate.
%  [ei,Tc] = extremalidxrange([Ie,xn(Ie,2)],opt);
%  subplot(2,1,1),plot(ei)
%  subplot(2,1,2),plot(Tc)
%  Tc = ceil(interp1(Tc.args,Tc.dataCI(:,2),1.2));
%  Tmin = xn(Tc,1)-xn(1,1);
%  Ie12 = findpot(xn,1.2,Tmin);
%  [ei12,Tc12] = extremalidxrange([Ie12,xn(Ie12,2)],opt,'umin',1.2);
%
% See also reslife, fitgenparrange, disprsnidxrange, findpot, decluster


% Reference
% Christopher A. T. Ferro, Johan Segers (2003)
% Inference for clusters of extreme values
% Journal of the Royal Statistical Society: Series B (Statistical Methodology) 65 (2), 545�556.
% doi:10.1111/1467-9868.00401

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


options = struct('u',[],'umin',[],'umax',[],'Nu',50,'Nmin',5,'Tb',1,'alpha',0.05,'plotflag',0,'method',6,'B',[]);

if ischar(data) && strcmpi(data,'defaults')
  ei = options;
  return
end
options = parseoptions(options,varargin{:});

[n m]  = size(data);
if min(m,n)==1,  % only data are given
  d  = data(:);
  n = length(d);
  t = 1:n;
else
  [t,ix] = sort(data(:,1)-min(data(:,1)));
  d = data(ix,2);
end


if isempty(options.u)
  sd = sort(d);
  if ~isempty(options.Nmin)
    options.Nmin = max(options.Nmin,0);
  end
  if options.Nmin>n/2
    warning('WAFO:EXTREMALIDXRANGE','Nmin possibly too large!')
  end
  if (n<=options.Nmin)
    error('Not enough data for a Extremal Index')
  end
  
  if isempty(options.umin)
    options.umin = sd(1);
  else
    options.umin = max(options.umin,sd(1));
   end
  
  sdmax = sd(n-options.Nmin);
  if isempty(options.umax)
    options.umax = sdmax;
  else
    options.umax = min(options.umax,sdmax);
  end

  if isempty(options.Nu)
    options.Nu = min(n-options.Nmin,100);
  end

  u = linspace(options.umin,options.umax,options.Nu).';
else
  u = options.u(:);
end

Nu = length(u);
nan1 = nan;
ei1 = nan1(ones(Nu,3));
TNmin = ei1;

for ix = 1:Nu
  excess = d > u(ix);
  [ei1(ix,:),TNmin(ix,:)] = ciextremalidx(t(excess),options.method,options.alpha,options.B);
end
p = 1-options.alpha;
eiCI = ei1(:,[1,3]);
tCI = TNmin(:,[1,3]);
  
CItxt = sprintf('%d%s CI',100*p,'%');
titleTxt = 'Extremal Index plot';

ei = createwdata('data',ei1(:,2),'args',u,...
'dataCI',eiCI,'title',titleTxt,'labels',{'Threshold','Extremal Index'},...
    'workspace',options,'note',titleTxt,'caption',CItxt);

Tc =  createwdata('data',TNmin(:,2),'args',u,...
'dataCI',tCI,'title','Minimum Cluster distance','labels',{'Threshold','distance'},...
    'workspace',options);
  

  ei = wdata(ei);
  Tc = wdata(Tc);
  if options.plotflag>0
    plot(ei,options.plotflag,'.')
  end


%%
function [ei,Tc] = ciextremalidx(t,method,alpha,B)
  Ti = diff(t); % interexceedance times
  ei = ciboot(Ti,@lcl_extrml_idx,method,alpha,B);

  if nargout>1
    N = length(Ti)+1;
    C = floor(N*ei)+1;
    sTi = -sort(-Ti);
    Tc = sTi(min(C,N-1)).'; % declustering time
    k = find(ei==1);
    if any(k)
      Tc(k) = min(Ti);    
    end
    Tc = fliplr(Tc);
  end


function ei = lcl_extrml_idx(Ti)
  Tmax = max(Ti); 
  if Tmax<=1,
    ei = 0;
  elseif Tmax<=2
    ei = min(1,2*mean(Ti).^2/mean(Ti.^2));
  else
    ei = min(1,2*mean(Ti-1).^2/mean((Ti-1).*(Ti-2)));
  end

