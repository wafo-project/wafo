function h = plotdensity(self)
%PLOTDENSITY Density plot.
%
% CALL:  H = densityplot(R,dist);
%
%   H    = handle to the plotted lines.
%   R    = data
%   dist = string containing the name of the PDF.
% 
%   DENSITYPLOT displays a density plot. The purpose of the plot is to 
%   graphically assess whether the data in R could come from the given 
%   distribution. If so the histogram should resemble the model density. 
%   Other distribution types will introduce deviations in the plot.  
%
%   This works on any PDF having the following calling syntax:
% 
%    phat = XXXfit(R);
%    p    = XXXpdf(p,phat(1),phat(2),...,phat(n));
%
%   where R contain the data and phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgenpar(1,2,0,100,1);
%   phat = fitgenpar(R);
%   phat.plotdensity; % or
%   plotdensity(phat);
%
% See also kdeplot, plotresp, plotresq, wqqplot

% Tested on: Matlab 5.3
% History:
%   Per A. Brodtkorb 12.11.2000


h = plotdensity(struct(self));


