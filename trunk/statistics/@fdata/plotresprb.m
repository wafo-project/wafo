function [varargout] = plotresprb(self)
%PLOTRESPRB Residual Probability plot.
%
% CALL:  H = plotresprb(self);
%
%   H    = handle to the plotted lines.
%   self = FDATA object with fields
%     .data         = data vector
%     .distribution = string containing the name of the PDF or function handle.
%     .params       = distribution parameters (default estimated from the R-data)
% 
%   PLOTRESPRB displays a residual probability plot. The purpose of the plot
%   is to graphically assess whether the data could come from the 
%   given distribution. If so the plot will be linear. Other distribution 
%   types will introduce curvature in the plot.  
%
%   This works on any PDF having the following calling syntax:
%    
%    p    = prbXX(R,phat(1),phat(2),...,phat(n))
%    phat = fitXX(R);
%
%   where R contain the data and phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgam(1,2,100,1);
%   plotresprb(R,'pdfgam');
%
% See also plotresq, plotqq

% Tested on: Matlab 5.3
% History:
% pab 2007
% renamed from distplot to resqplot
% by  Per A. Brodtkorb 12.11.2000

error(nargchk(1,1,nargin))
[varargout{1:nargout}]=plotresprb(struct(self));
return
% dist = self.distribution;
% x    = self.data{1};
% phat = self.params;
% if isa(dist,'funciton_handle')
%   dist = func2str(dist);
% end
% 
% if strncmpi(fliplr(dist),'fdp',3) || strncmpi(fliplr(dist),'vni',3) || strncmpi(fliplr(dist),'tif',3)
%   pdf = dist(1:end-3);
% else
%   pdf = dist; 
% end
% 
% x=sort(x(:));
% 
% 
% if isempty(phat)
%   self  = feval( [ pdf 'fit'],x);% MLE of the distribution parameters
%   phat = self.params;
% end
% 
% cphat = num2cell(phat,1);
% 
% n     = length(x);
% ecdf = (0.5:n-0.5)/n;
% mcdf  = feval([ pdf 'cdf'],x,cphat{:})';
% 
% if nargout > 0
%    h = wqqplot(ecdf,mcdf);
% else 
%    wqqplot(ecdf,mcdf);
% end
% 
% xlabel('Empirical')
% ylabel(sprintf('Model: %s',pdf))
% title('Residual Probability Plot');
% axis('equal')
% wafostamp;
% %grid on;
% 
