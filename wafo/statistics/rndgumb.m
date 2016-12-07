function R  = rndgumb(varargin)
%RNDGUMB Random matrices from a Gumbel distribution.
% 
% CALL:  R = rndgumb(a,b,sz,options) 
%  
%  R     = a matrix of random numbers from the Gumbel distribution
%  a, b  = parameters of the Gumbel distribution.
%   phat = Distribution parameter struct
%             as returned from FITGUMB.  
%     sz = size(R)    (Default common size of k,s and m0)
%          sz can be a comma separated list or a vector 
%          giving the size of R (see zeros for options). 
%   options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .trunc    : if TRUE truncated gumbel distribution
%                       otherwise regular gumbel distribution (default)
%
%   The size of R is the common size of a and b if both are matrices.
%   If either parameter is a scalar, the size of R is the size of the other
%   parameter. R = rndgumb(a,b,trunc,m,n) returns an m by n matrix. 
%
% Example:
%   R=rndgumb(5,10,1,100);
%   plotgumb(R);
%
%   close all;
%
% See also  fitgumb, pdfgumb, cdfgumb, invgumb, momgumb

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


% Reference: 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 

% Tested on: matlab 5.2
% History:
% revised pab 8.11.1999
%  updated help header
% by Per A. Brodtkorb 26.10.98
% revised ms 13.06.2000
% - updated header info
% - changed name to rndgumb (from gumbrnd)
% - added w* to used WAFO-files 
% - enabled use of 5:th argument
% - removed stat-toolbox routines
% 
% revised pab 23.10.2000
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R

error(nargchk(1,inf,nargin))
Np = 2;
options = struct('lowertail',true,'trunc',false); % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
b = params{2};

if isempty(options.trunc)
  options.trunc=false; %default is not truncated gumbel
end
if isempty(rndsize)
  csize = comnsize(a,b);
else
  csize = comnsize(a,b,zeros(rndsize{:}));
end
if any(isnan(csize))
  error('a and b must be of common size or scalar.');
end
R=invgumb(rand(csize),a,b,options);




