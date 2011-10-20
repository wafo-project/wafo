function h = plotresq(x,dist, varargin)
%PLOTRESQ Plot Residual Quantile.
%
% CALL:  H = plotresq(R,dist,phat(1),phat(2),...,phat(n));
%        H = plotresq(phat1);
%
%   H    = handle to the plotted lines.
%   R    = data
%   dist = string containing the name of the PDF or function handle.
%  phat1 = fdata struct of distribution parameters and data.
%   phat = vector of distribution parameters (default estimated from the R-data)
% 
%   PLOTRESQ displays a residual quantile plot. The purpose of the plot
%   is to graphically assess whether the data in R could come from the 
%   given distribution. If so the plot will be linear. Other distribution 
%   types will introduce curvature in the plot.  
%
%   This works on any PDF having the following calling syntax:
%    
%    p    = XXXpdf(R,phat(1),phat(2),...,phat(n))
%    phat = XXXfit(R);
%    x    = XXXinv(p,phat(1),phat(2),...,phat(n));
%
%   where R contain the data and phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgam(1,2,100,1);
%   plotresq(R,'pdfgam');
%
% See also plotresprb plotqq

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


% Tested on: Matlab 5.3
% History:
% pab 2007
% renamed from distplot to resqplot
% by  Per A. Brodtkorb 12.11.2000

error(nargchk(1,inf,nargin))
cphat = [];
if isa(x,'struct') || isa(x,'fdata')
  ih = ishold;
  Np = numel(x);
  if Np>1
    cfig=gcf;
    H1 = cell(1,Np);
    for ix = 1:Np
       if ih
        newplot
      else
        figure(cfig-1+ix)
      end
      H1{ix} = plotresq(x(ix));
    end
    if nargout>0
      h = H1;
    end
    return
  end
  cphat  = num2cell(x.params,1);
  dist  = x.distribution;
  x = x.data(:);
else
  error(nargchk(2,inf,nargin))
  x=x(:); 
  if nargin>2
    cphat = varargin;
  end
end
  
model = getdistname(dist);
if isempty(cphat)
  phat  = feval( ['fit' model],x);% MLE of the distribution parameters
  cphat = num2cell(phat,1);
end

n=length(x);
eprob = ((1:n)-0.5)/n;
y  = feval(['inv' model],eprob,cphat{:})';

if nargout > 0
   h = plotqq(x,y);
else 
   plotqq(x,y);
end

xlabel('Empirical')
ylabel(sprintf('Model (%s)',model))
title('Residual Quantile Plot');
axis('equal')
axis('square')
wafostamp;

if 0 
  % Old call kept just in case
p     = [0.001 0.02 0.05 0.10 0.25 0.5 0.75 0.90 0.95 0.98 0.99 0.997 0.999];
tick  = feval(['inv' model ],p,cphat{:});
ax=axis;hold on
plot([ax(1) ax(2)],[tick; tick],'b:'); hold off
for l=1:length(p)
  h1=figtext(1.03,tick(l),num2str(p(l)) ,'norm','data');
  set(h1,'FontSize',10,'FontAngle','Italic')
end
end
%grid on;

