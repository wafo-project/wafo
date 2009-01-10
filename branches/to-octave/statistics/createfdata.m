function f=createfdata(varargin)
% CREATEFDATA Distribution parameter struct constructor
%
% CALL:  f=createfdata(list)
%
%  f = struct with the following fields
%   .distribution : function handle or name of PDF function
%   .params       : vector of distribution parameters
%   .lowerbound   : Lower bound of PARAMS estimate at 100*(1-alpha)% level
%   .upperbound   : Upper bound of PARAMS estimate at 100*(1-alpha)% level
%   .fixpar        : vector giving vector giving the fixed parameters.
%   .alpha        : Confidence coefficient
%   .covariance   : Covariance matrix of estimate
%   .variance     : diag(covariance)
%   .loglikemax   : maximum of log-likelihood
%   .logpsmax     : maximum of Moran's log product spacing statistic
%   .pvalue       : p-value for the hypothesis, H0, that model fits data.
%   .dataname     : dataname
%   .data         : cellarray with data if params is estimated
%   .method       : method used in the parameter estimation
%   .pdfoptions   : options to submit to pdfXX, prbXX
%   .note         : note string
%   .date         : creation date and time
%
%
%  Examples:   
%   f = createfdata('dist','pdfgenpar')  %gives the structure
% %     distribution: 'pdfgenpar'
% %           params: []
% %       lowerbound: []
% %       upperbound: []
% %            alpha: []
% %       covariance: []
% %         variance: []
% %       loglikemax: []
% %             data: []
% %           method: []
% %             note: []
% %             date: '07-Sep-2007 11:09:24'
%
% See also  fdata

%Tested on: Matlab 5.3
%History:
% by pab Sep 2007



names = {'distribution','params','upperbound','lowerbound',...
  'fixpar','alpha', 'covariance','variance','loglikemax',...
  'logpsmax','pvalue','dataname','data','method','pdfoptions','note','date'};

% names = {'data','distribution','method',...
%   'phat','phat_lower','phat_upper','alpha',...
%   'pcov','pvar',...
%   'loglike_max','note','date'};

n    = length(names);
c    = cell(1,n);
c{n} = datestr(now);
f    = cell2struct(c,names,2);
if nargin>0
  f = parseoptions(f,varargin{:});
end

