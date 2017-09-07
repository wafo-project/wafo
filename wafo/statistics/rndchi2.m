function R = rndchi2(varargin)
%RNDCHI2 Random matrices from a Chi squared distribution.
%
%  CALL:  R = rndchi2(df,sz);
%
%        df = degrees of freedom, df=1,2,3,...
%      phat = Distribution parameter struct
%             as returned from FITCHI2.  
%        sz = size(R)    (Default size(df))
%             sz is a comma separated list or a vector 
%             giving the size of R (see zeros for options).
%
%  The Chi squared distribution is a special case of 
%  the gamma distribution. Thus the rndgam is used to
%  generate R.
%
% Example:
%   R=rndchi2(2,1,100);
%   plot(R,'.');
%
%   close all;
%
% See also  rndgam, pdfchi2, cdfchi2, fitchi2, momchi2

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
% revised pab Dec2003
% removed old call
% revised pab 25.10.2000
% - replaced code with a call to rndgam
% added ms 26.06.2000

%error(nargchk(1,inf,nargin))
narginchk(1,inf)
Np = 1;
options = []; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
% if numel(options)>1
%   error('Multidimensional struct of distribution parameter not allowed!')
% end

df = params{1};
ok = (df==round(df)& 0<df);
if any(~ok)
  warning('WAFO:RNDCHI2','df must be a positive integer')
  df(~ok) = nan;
end
if isempty(rndsize),
  csize = size(df);
else
  [csize,df] = comnsize(df,zeros(rndsize{:}));
end

if any(isnan(csize))
  error('df must be a scalar or comply to the size given.');
end

R = rndgam(df/2,2);





