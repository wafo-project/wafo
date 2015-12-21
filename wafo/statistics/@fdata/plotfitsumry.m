function H = plotfitsumry(self,varargin)
% PLOTFITSUMRY Diagnostic plot for the fit to data
% 
% CALL h = plotfitsummary(self,plotflag)
%
% H    = handle to the plotted lines.
% data = vector
% dist = string containing the name of the PDF or function handle.
% phat = distribution parameters (default estimated from the data)
% plotflag = see empdistr for details.
%  
%    FITSUMRYPLOT displays empirical CDF and PDF vs model as well as 
%    kernel density plot, residual quantile plot and residual probability 
%    plot. The purpose of the plots is to graphically assess whether the 
%    data in R could come from the given distribution. If so the empirical- 
%    CDF and PDF should follow the model and the residual plots will be linear.
%    Other distribution types will introduce curvature in the residual plots.  
%
% Example
%  R = rndgenpar(1,2,0,100,1);
%  phat = fitgenpar(R);
%  phat.fitsumryplot; % or
%  fitsumryplot(phat);
%
% See also edf, plotresq, plotresp, plotqq


H = plotfitsumry(struct(self),varargin{:});
