function [param,options,rndsize,phat] = parsestatsinput(Np,options,varargin)
%PARSESTATSINPUT Parses inputs to pdfxx, prbxx, invxx and rndxx functions.
%
% CALL  [param, options,sz] = parsestatsinput(N,options0,a1,a2,...,an,size,options1)
%       [param, options,sz] = parsestatsinput(N,options0,phat,size,options1)
%       [param, options,sz] = parsestatsinput(N,options0,phat,size,options1)
%
%  param      = cellarray of distribution parameters size 1 x N
%  options    = merged options structure
%  sz         = cellarray with size info
%  N          = Number of parameters
%  options0   = default options struct.
%  a1,a2,.,an = Distribution parameter values
%  phat       = Distribution parameter struct
%  size       = vector or list with sizeinfo.
%  options1   = options structs or lists of named parameters and
%               corresponding values.
%
% Example
%  opt = struct('covariance',[],'alpha',[]);
%  [param,options] = parsestatsinput(2,opt,1,'cov',2,'alpha',0.05);
%  assert(options, struct('covariance', 2, 'alpha', 0.05))
%  assert(param, {1,[]})
%
% See also parseoptions


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


%error(nargchk(2,inf,nargin))
narginchk(2,inf)
ni = nargin-2;
if isnan(Np), Np = 30;end
param   = cell(1,Np);
rndsize = cell(1,0);
phat = [];
if ni>0
  isnum = cellfun(@isnumeric,varargin);
  
  N = find(~isnum,1,'first');
  if isempty(N)
    N = length(isnum)+1;
  end
  if (N==1 && (isa(varargin{1},'struct') && all(isfield(varargin{1},{'params','distribution'})) || isa(varargin{1},'fdata')))  
    try
      phat = struct(varargin{1}); % fit struct
      param1 = num2cell(cat(1,phat.params),1);
      Np1 = min(length(param1),Np);
      param(1:Np1) = param1(1:Np1);
      if nargout>1
        nistart = find(~isnum(2:ni),1,'first')+1;
        if isempty(nistart)
          % xxxrnd(phat,size)
          % xxxpdf(phat)
          nistart = length(isnum)+1;
          % else
          % xxxrnd(phat,size,options)
          % xxxpdf(phat,options)
        end
      
        rndsize = varargin(2:nistart-1);
        if ~isempty(options)
          nump =  numel(phat);
          if nump>1 % expand
            [options(1:nump) ] = deal(options);
          end
          if isfield(phat,'pdfoptions') && ~isempty(phat.pdfoptions)
            options = parseoptions(options,phat.pdfoptions,struct('alpha',phat.alpha),varargin{nistart:ni});
          else
            options = parseoptions(options,struct('alpha',phat.alpha),varargin{nistart:ni});
          end
        end
      end
    catch
      error('WAFO:PARSESTATSINPUT','First input must either be a numeric value or a distribution parameter struct.')
    end
  else
    % xxxpdf(a1,a2,...,an1,options)
    % xxxrnd(a1,a2,...,anp,rndsize,options)
    N1 = min(Np,N-1);
    param(1:N1)       = varargin(1:N1);
    if nargout>1
      rndsize           = varargin(N1+1:N-1);
      if ~isempty(options)
        options           = parseoptions(options,varargin{N:ni});
      end
    end
  end
end

