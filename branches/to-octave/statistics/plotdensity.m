function h = plotdensity(x,dist,varargin)
%PLOTDENSITY Plot density.
%
% CALL:  H = plotdensity(R,dist,phat(1),phat(2),...,phat(n));
%        H = plotdensity(phat1);
%
%   H    = handle to the plotted lines.
%   R    = data
%   dist = string containing the name of the PDF or function handle.
%  phat1 = fdata struct of distribution parameters and data.
%   phat = distribution parameters (default estimated from the R-data)
% 
%   PLOTDENSITY displays a density plot. The purpose of the plot is to 
%   graphically assess whether the data in R could come from the given 
%   distribution. If so the histogram should resemble the model density. 
%   Other distribution types will introduce deviations in the plot.  
%
%   This works on any PDF having the following calling syntax:
% 
%    phat = fitXXX(R);
%    p    = pdfXXX(p,phat(1),phat(2),...,phat(n));
%
%   where R contain the data and phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgam(1,2,100,1);
%   plotdensity(R,'pdfgam');
%
% See also plotkde, plotresprb, plotresq, plotqq

% Tested on: Matlab 5.3
% History:
%   Per A. Brodtkorb 12.11.2000

ih = ishold;
cphat = [];
if isa(x,'struct') || isa(x,'fdata')
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



[r,xmax,xmin] = range(x);
rmin = xmin-r/4;
rmax = xmax + r/4;
fhandle = str2func(['pdf' model]);
f1 = @(x) fhandle(x,cphat{:});
[x1,y1] =  fplot(f1,[rmin,rmax]);


histgrm(x,[],0,1)
hold on
h1 = plot(x1,y1,'b-');
hold off

if nargout > 0
   h = h1;
end


xlabel('x');
ylabel(sprintf('f(x) (%s)',model)) 
title('Density plot');

wafostamp;


