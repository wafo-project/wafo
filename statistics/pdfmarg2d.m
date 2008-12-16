function f = pdfmarg2d(V,H,varargin)
%PDFMARG2D Joint 2D PDF due to Plackett given as  f{x1}*f{x2}*G(x1,x2;Psi). 
%
%  CALL:  f = pdfmarg2d(x1,x2,phat,options) 
%
%     f  = PDF evalutated at points (x1 , x2) 
%   phat = parameter structure as returned from fitmarg2d
%         .params : [par1,par2,psi] Vector where 
%                    par1 is parameters for X1         
%                    par is parameters for X2
%                    psi is interaction parameter between x1 and x2.
% options = options structure with fieldnames:
%    .distribution : 2D cellarray of distribution names
%    .numpar = vector of number of parameters
%    .condon = 0 regular pdf is returned (default)
%              1 conditional pdf of X2 given X1 is returned
%              2 conditional pdf of X1 given X2 is returned
%    .meshgrid : if true compute f on meshgrid(x1,x2) (default false)
%    .wdata    : if true return f as a wdata object   (default false)
%
% Example: % 2D Weibull Rayleigh with parameters [2 3] and 3,
%          % respectively for the marginal distributions, and
%          % a interaction parameter of 10:
%  opts = pdfmarg2d('defaults');
%  params = {2,3,3,10}
%  opts=parseoptions(opts,'numpar',[2,1],'distribution',{'pdfweib','pdfray'},'meshgrid',true,'wdata',true);
%  x = linspace(0,5,50); x2 = linspace(0,10);
%  f = pdfmarg2d(x,x2,params{:},opts);
%  plot(f)
% 
% See also   cdfmarg2d, invcmarg2d, fitmarg2d, rndmarg2d


%   References:
%      Plackett, R. L. (1965) "A class of bivariate distributions."
%                                J. Am. Stat. Assoc. 60. 516-22
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 pp. 133-134.


%  tested on: matlab 5.2
% history
% revised pab 20.10.2000
% - updated to new wstats toolbox
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vectro to structure
%  Per A. Brodtkorb 28.01.99


options = struct('numpar',[],'distribution',[],'condon',0,...
  'meshgrid',false,'wdata',false,'logp',false);
if nargin==1 && nargout <= 1 && isequal(V,'defaults')
  f = options;
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
  error ('x1 and x2 must be of common size or scalar');
end
try
  VDIST=getdistname(lower(options.distribution{1}));
  HDIST=getdistname(lower(options.distribution{2}));
catch
  if isempty(phat)
    error('Not enough inputs')
  end
  VDIST=getdistname(lower(phat.pdfoptions.distribution{1}));
  HDIST=getdistname(lower(phat.pdfoptions.distribution{2}));
end

% psi = options.psi; % interaction parameter
% PV  = num2cell(options.param1(:).',1);
% PH  = num2cell(options.param2(:).',1);


Fvh = str2func(['cdf',VDIST]);
fvh = str2func(['pdf',VDIST]);
%mvh = str2func(['mom',VDIST]);
Fhh = str2func(['cdf',HDIST]);
fhh = str2func(['pdf',HDIST]);
%mhh = str2func(['mom',HDIST]);
  
Fv = Fvh(V,PV{:});
fv = fvh(V,PV{:});
Fh = Fhh(H,PH{:});
fh = fhh(H,PH{:});
%  Fv=dist1dcdffun(V(k),PV, VDIST(1:2));
%  fv=dist1dpdffun(V(k),PV, VDIST(1:2));
%  Fh=dist1dcdffun(H(k),PH, HDIST(1:2));
%  fh=dist1dpdffun(H(k),PH, HDIST(1:2));
 
  tmp=1+(Fv+Fh).*(psi-1);
  f = psi.*((psi-1).*(Fv+Fh-2.*Fv.*Fh)+1)./(sqrt(tmp.^2-4.*psi.*(psi-1).*Fv.*Fh).^3);
  k = find(isfinite(f) & f~= 0);
  if any(k)
   switch options.condon
     case 0, f(k)=f(k).*fv(k).*fh(k);
     case 1, f(k)=f(k).*fh(k);
     case 2, f(k)=f(k).*fv(k);
     case 3, % secret option  used by mommarg2d: returns v*f(v|h)
       f(k)=V(k).*f(k).*fv(k);
     case 4, % secret option  used by mommarg2d: returns v^2*f(v|h)
       f(k)=V(k).^2.*f(k).*fv(k);
     case 13, % secret option  used by mommarg2d: returns h*f(h|v)
       f(k)=H(k).*f(k).*fh(k);
     case 14, % secret option  used by mommarg2d: returns h^2*f(h|v)
       f(k)=H(k).^2.*f(k).*fh(k);
   end
  end
  if options.logp
    f = log(f);
  end
%plot(V(k),y(k),'.') 
 if options.wdata
   try
     [cl pl]=qlevels(f);
   catch
     cl = [];
   end
   captn = sprintf('pdfmarg2d: (%s,%s)',VDIST,HDIST);
   f = createwdata('args',{Vin,Hin},'data',f,'caption',captn,'contourLevels',cl,'percentLevels',pl,'workspace',options);
   if ~isoctave
     f = wdata(f);
   end
 end
end







