function [varargout] = plotresq(self)
%PLOTRESQ Residual Quantile plot.
%
% CALL:  H = plotresq(self);
%
%   H    = handle to the plotted lines.
%   self = FDATA object with fields
%     .data         = data vector
%     .distribution = string containing the name of the PDF or function handle.
%     .params       = distribution parameters (default estimated from the R-data)
% 
%   PLOTRESQ displays a residual quantile plot. The purpose of the plot
%   is to graphically assess whether the data could come from the 
%   given distribution. If so the plot will be linear. Other distribution 
%   types will introduce curvature in the plot.  
%
%   This works on any PDF having the following calling syntax:
%    
%    p    = pdfXX(x,phat(1),phat(2),...,phat(n))
%    phat = fitXX(data);
%    x    = invXX(p,phat(1),phat(2),...,phat(n));
%
%   where phat(1),phat(2)... are
%   the distribution parameters. 
%
% Example:
%   R = rndgam(1,2,100,1);
%   plotresq(R,'pdfgam');
%
% See also plotresprb plotqq

% Tested on: Matlab 5.3
% History:
% pab 2007
% renamed from distplot to resqplot
% by  Per A. Brodtkorb 12.11.2000

%error(nargchk(1,1,nargin))
narginchk(1,1)
[varargout{1:nargout}]=plotresq(struct(self));
return


% dist = self.distribution;
% if iscell(self.data)
%   x = self.data{1};
% else
%   x = self.data;
% end
% phat = self.params;
% if isa(dist,'function_handle')
%   dist = func2str(dist);
% end
% 
% model = getdist
% if strncmpi(fliplr(dist),'fdp',3) || strncmpi(fliplr(dist),'vni',3) || strncmpi(fliplr(dist),'tif',3)
%   pdf = dist(1:end-3);
% else
%   pdf = dist; 
% end
% 
% x=x(:);
% 
% 
% if isempty(phat)
%   self  = feval( [ pdf 'fit'],x);% MLE of the distribution parameters
%   phat = self.params;
% end
% 
% cphat = num2cell(phat,1);
% 
% n=length(x);
% eprob = ((1:n)-0.5)/n;
% y  = feval([ pdf 'inv'],eprob,cphat{:})';
% 
% if nargout > 0
%    h = wqqplot(x,y);
% else 
%    wqqplot(x,y);
% end
% 
% xlabel('Empirical')
% ylabel(sprintf('Model: %s',pdf))
% title('Residual Quantile Plot');
% axis('equal')
% wafostamp;
% 
% %grid on;
% 
