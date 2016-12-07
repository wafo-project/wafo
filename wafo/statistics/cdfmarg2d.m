function F = cdfmarg2d(V,H,varargin)
%CDFMARG2D Joint 2D CDF due to Plackett 
%
% CALL F = cdfmarg2d(x1,x2,phat,options) 
%
%     F  = the cdf evalutated at points (x1 , x2) 
%          with the parameters Phat1, Phat2 and Psi.
%   phat = parameter structure as returned from fitmarg2d.
%         .params : [param1,param2,psi] where param1 and param2 are vectors of
%                    marginal parameters for X1 and X2 respectively. 
%                    psi is the interaction parameter between x1 and x2.
% options = options structure with fieldnames:
%    .distribution : 2D cellarray of distribution names
%    .numpar = vector of number of parameters
%    .condon = 0 regular cdf is returned (default)
%              1 conditional cdf of X2 given X1 is returned
%              2 conditional cdf of X1 given X2 is returned
%    .meshgrid : if true compute f on meshgrid(x1,x2) (default false)
%    .wdata    : if true return f as a wdata object   (default false)
%
% Example: % 2D Weibull Rayleigh with marginal parameters [2 3] and 3,
%         % respectively and interaction parameter of 10 : 
%  opts = cdfmarg2d('defaults');
%  params = {2,3,3,10};
%  opts = parseoptions(opts,'numpar',[2,1],'distribution',...
%                      {'pdfweib','pdfray'},'meshgrid',true,'wdata',true);
%  x = linspace(0,5,50); x2 = linspace(0,10);
%  F = cdfmarg2d(x,x2,params{:},opts);
%  surf(F);
%
%  close all;
% 
% See also pdfmarg2d, fitmarg2d, rndmarg2d

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



% references:
% Plackett, R. L. (1965) "A class of bivariate distributions."
%                                J. Am. Stat. Assoc. 60. 516-22
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 p. 133-134.

%  tested on: matlab 5.2
% history
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vectro to structure
%  Per A. Brodtkorb 28.01.99

% options = struct('param1',[],'param2',[],'psi',[],'dist1','','dist2','',...
%   'condon',0, 'meshgrid',false,'wdata',false);
options = struct('numpar',[],'distribution',[],'condon',0,...
  'meshgrid',false,'wdata',false);

if nargin==1 && nargout <= 1 && isequal(V,'defaults')
  F = options;
  return
end
error(nargchk(3,inf,nargin))
Np = nan;

[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
PV = params(1:options.numpar(1));
PH = params(options.numpar(1)+(1:options.numpar(2)));
psi = params{sum(options.numpar)+1};


Vin = V;
Hin = H;
if options.meshgrid
  [V,H] = meshgrid(V,H);
end


[icode V H ] = iscomnsize(V,H);
if  ~icode
  error('x1 and x2 must be of common size or scalar.');
end

try
  VDIST=getdistname(lower(options.distribution{1}));
  HDIST=getdistname(lower(options.distribution{2}));
catch
  if isempty(phat)
    error('Not enough inputs')
  end
  VDIST=getdistname(lower(phat.distribution{1}));
  HDIST=getdistname(lower(phat.distribution{2}));
end


% VDIST=getdistname(lower(options.dist1));
% HDIST=getdistname(lower(options.dist2));


% psi = options.psi; % interaction parameter
% PV  = num2cell(options.param1(:).',1);
% PH  = num2cell(options.param2(:).',1);


Fvh = str2func(['cdf',VDIST]);
Fhh = str2func(['cdf',HDIST]);


Fv = Fvh(V,PV{:});
Fh = Fhh(H,PH{:});

tmp=1+(Fv+Fh).*(psi-1);
switch options.condon
  case 0, F = (tmp-sqrt(tmp.^2-4.*psi.*(psi-1).*Fv.*Fh))./(2.*(psi-1));
  case 1, F = 0.5-0.5.*(tmp-2.*psi.*Fh)./sqrt(tmp.^2-4.*psi.*(psi-1).*Fv.*Fh);
  case 2, F = 0.5-0.5.*(tmp-2.*psi.*Fv)./sqrt(tmp.^2-4.*psi.*(psi-1).*Fv.*Fh);
end


if options.wdata
  captn = sprintf('cdfmarg2d: (%s,%s)',VDIST,HDIST);
  cl = [];
  pl = [];
  F = createwdata('args',{Vin,Hin},'data',F,'caption',captn,'contourLevels',cl,'percentLevels',pl,'workspace',options);
  if ~isoctave
    F = wdata(F);
  end
end

end


