function [x,xlo,xup] = invchi2(F,varargin)
%INVCHI2 Inverse of the Chi squared distribution function
%
% CALL:  x = invchi2(F,df,options)
%        [x,xlo,xup] = invchi2(F,phat,options);
%
%        x = inverse cdf for the Chi squared distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%       df = degrees of freedom
%     phat = Distribution parameter struct
%            as returned from FITCHI2.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 500)
%
%
% Example:
%     F = linspace(0,1,100);
%     x = invchi2(F,1);
%     plot(F,x);
%
%     close all;
%
% See also invgam, pdfchi2, cdfchi2, rndchi2, fitchi2, momchi2

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


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 415 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% revised pab 25.10.2000
%  - added comnsize, nargchk
% added ms 26.06.2000
%error(nargchk(2,inf,nargin))
narginchk(2,inf)
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',500); % default options

Np = 1;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

df = params{1};
if any(df~=round(df)),
  warning('WAFO:INVCHI2','df should be a positive integer')
end
try
  if nargout>1
    if isempty(phat)
      error('Must have distribution struct!')
    end
    phat.covariance = [phat.covariance/4,0:0 0];
    phat.distribution = 'invf';
    phat.params = [df/2,2];
    phat.fixpar = [phat.fixpar/2,2];
    [x,xlo, xup] = invgam(F,phat,options);
  else
    x=invgam(F,df/2,2,options);
  end
catch
  error('F and df must be of common size or scalar.');
end


