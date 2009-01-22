function h = plotresprb(x,dist,varargin)
%PLOTRESPRB Plot Residual Probability.
%
% CALL:  H = plotresprb(R,dist,phat(1),phat(2),...,phat(n));
%        H = plotresprb(phat1);
%
%   H    = handle to the plotted lines.
%   R    = data
%   dist = string containing the name of the PDF or function handle.
%  phat1 = fdata struct of distribution parameters and data.
%   phat = distribution parameters (default estimated from the R-data)
% 
%   PLOTRESPRB displays a residual probability plot. The purpose of the plot
%   is to graphically assess whether the data in R could come from the 
%   given distribution. If so the plot will be linear. Other distribution 
%   types will introduce curvature in the plot.  
%
%   This works on any PDF having the following calling syntax:
%    
%    p    = XXXcdf(R,phat(1),phat(2),...,phat(n))
%    phat = XXXfit(R);
%
%   where R contain the data and phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgam(1,2,100,1);
%   plotresprb(R,'pdfgenpar');
%
% See also plotresq plotqq

% Tested on: Matlab 5.3
% History:
% pab 2007
% renamed from distplot to resqplot
% by  Per A. Brodtkorb 12.11.2000

error(nargchk(1,inf,nargin))


cphat = [];
if isa(x,'struct') || isa(x,'fdata')
  Np = numel(x);
  if Np>1
    ih = ishold;
    cfig = gcf;
    H1 = cell(1,Np);
    for ix = 1:Np
      if ih
        newplot
      else
        figure(cfig-1+ix)
      end
      H1{ix} = plotresprb(x(ix));
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
  cphat = num2cell(phat.params,1);
end

n     = length(x);
%ecdf = (0.5:n-0.5)/n;
ecdf = (1:n)/(n+1);
mcdf  = feval(['cdf' model],x,cphat{:})';

if nargout > 0
   h = plotqq(ecdf,mcdf);
else 
   plotqq(ecdf,mcdf);
end

xlabel('Empirical')
ylabel(sprintf('Model (%s)',model))
title('Residual Probability Plot');
axis([0 1 0 1])
axis('square')
wafostamp;
%grid on;

