function H = plotfitsumry(phat,plotflag)
%PLOTFITSUMRY Plot diagnostic of fit to data
% 
% CALL h = plotfitsumry(phat,plotflag)
%
% H    = handle to the plotted lines.
% data = vector
% dist = string containing the name of the PDF or function handle.
% phat = distribution parameters (default estimated from the data)
% plotflag = see edf for details.
%  
%    PLOTFITSUMRY displays probability plot, density plot, residual quantile
%    plot and residual probability plot. The purpose of these plots is to 
%    graphically assess whether the data in R could come from the given 
%    distribution. If so the empirical- CDF and PDF should follow the model
%    and the residual plots will be linear. Other distribution types will 
%    introduce curvature in the residual plots.  
%
% Example
%  R = rndgam(1,2,100,1);
%  phat = fitgam(R);
%  plotfitsumry(phat);
%
% See also edf, plotresq, plotresprb, plotqq

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


%error(nargchk(1,2,nargin))
narginchk(1,2)
if nargin<2 || isempty(plotflag)
  plotflag = 1011;
end
hold_state = ishold; % remember old hold state
Nff=length(phat);
if Nff>1
  cfig=gcf;
  H = cell(1,Nff);
  for ix=1:Nff,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    H{ix} = plotfitsumry(phat(ix),plotflag);
  end
  return
end

dist = phat.distribution;
model = getdistname(dist);

pdf  = ['pdf' model];
cdf  = ['cdf' model];
fit =  ['fit' model];
data=sort(phat.data(:));

if isempty(phat.params)
  phat  = feval(fit,data);% MLE of the distribution parameters
end
cphat = num2cell(phat.params,1);
x = data;

f0 = feval(pdf,x,phat);
try
  [F0,Flo,Fup] = feval(cdf,x,phat);%,'proflog',false);
catch
  [F0] = feval(cdf,x,phat);
  Flo = F0*nan;
  Fup = F0*nan;
end

subplot(2,2,1)
plotedf(data,[x, F0],plotflag)

if any(isfinite(Flo)) || any(isfinite(Fup))
  Fci = wdata(F0,x);
  Fci = set(Fci,'dataCI',[Flo(:) Fup(:)]);
  hold on
  plot(Fci,plotflag)
  hold off
end
title('Probability plot')


switch lower(pdf)
  case {'pdfgenpar','pdfchi2','pdfexp','pdff','pdfgam','pdfgengam','pdfinvnorm','pdfray','pdflognorm','pdfraymod','pdfweib','pdfweibmod'}
    L2 = 0.5;
  otherwise
    L2 = 1;% No transformation
end

subplot(2,2,2)
histgrm(x,[],.0,1)
ax = axis;
hold on
f = kdebin(x,{'L2',L2,'hsMethod','hste'});
h2 = plot(x,f0,'k-',f.x{:},f.f,'b:');
hold off
ax2 = axis;
ax2(1:2) = ax(1:2);
if ax2(4)>2*ax(4)
  ax2(4) = min(ax2(4),1.5*ax(4));
  axis(ax2)
else
  axis('tight')
end
title('Density plot')
xlabel('x')

subplot(2,2,3)
h3 = plotresq(data,dist,cphat{:});
%ax3 = axis;
%ix = [1 3];
% ax3(ix) = (min(ax3(ix)));
% ax3(ix+1) = (max(ax3(ix+1)));
%axis(ax3)
axis('tight')

subplot(2,2,4)
h4 = plotresprb(data,dist,cphat{:});
axis([0 1 0 1])
if nargout > 0
   H = [h2(:);h3(:);h4(:)];
end

if isempty(phat.method)
  methodstr = '';
else
  methodstr = phat.method;
end
if ~isempty(phat.fixpar)
  isfixed = isfinite(phat.fixpar);
  i_fixed  = find(isfixed);
  if any(i_fixed)
    phatistr = sprintf('%d,',i_fixed);
    phatistr(end) = '';
    phatvstr = sprintf('%g,',phat.params(i_fixed));
    phatvstr(end) = '';
    if numel(isfixed)>1
      phatvstr = sprintf('[%s]',phatvstr);
    end
    fixstr = sprintf('Fixed: phat(%s) = %s ',phatistr,phatvstr);
  else
    fixstr = '';
  end
else
  fixstr = '';
end

infostr = sprintf('Fit method: %s, Fit p-value: %2.2f %s',methodstr,phat.pvalue,fixstr);
 wafostamp(infostr)
