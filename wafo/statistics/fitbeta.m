function [phat]=fitbeta(data, varargin)
%FITBETA Parameter estimates for Beta data.
%
% CALL:  [phat] = fitbeta(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ml' Maximum Likelihood method (default)
%               'mps' Maximum product of spacings method
%               'mom' Moments method
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%          
% Example:
%   R=rndbeta(2,2,1,100);
%   phat = fitbeta(R)
%   plotfitsumry(phat)
%
% See also  pdfbeta, cdfbeta, rndbeta, invbeta, mombeta

% Copyright (C) 2000 Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
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

% No Reference
%  

% tested on: matlab 5.3
% History:
% By  pab 24.10.2000
% revised PJ 03-Apr-2001
%  - fmins changed name to fminsearch for version >= 5.3

global WAFO_WSTATS_DEFAULT_PLOTFLAG

error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

if any(data(:)<=0)
  error('data must be strictly positive!')
end


data=data(:); % make sure it is a vector

mu = mean(data);
sa = std(data)^2;

% Supply a starting guess with method of moments:
a = (mu*(1-mu)/sa-1)*mu;
phat0 = [ a a*(1/mu-1)];

ismom =  strcmpi(options.method,'mom');
if strcmpi(options.method,'ml') || strcmpi(options.method,'mps') || ismom
  % Maximum Likelihood , Maximum Product of Spacing or Moments estimator
  pdf = @pdfbeta;
  phat = mlest(pdf,phat0,data,{options,'search',~ismom});
  phat.dataname = inputname(1);
else
  error(['Unknown method ' options.method '.']);
end
