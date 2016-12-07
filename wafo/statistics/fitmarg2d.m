function  phato=fitmarg2d(V,H,varargin)
% FITMARG2D Parameter estimates for MARG2D data.
%
%  CALL: Phat = fitmarg2d(x1,x2,options)
%
%   phat = structure array containing estimated parameters
% x1,x2  = input data
% options = struct with fieldnames
%       .distribution : list of distributions used in the fitting of X1 and X2, 
%                       respectively.
%       .method : a string, describing the method of estimation
%                'sml'  = Simultainous Maximum Likelihood method (default)
%                'mml'  = Marginal Maximum Likelihood method (default)
%                'mmps' = Marginal Maximum Product Spacing method
%       .start   : startvalues used in fitting
%       .fixpar  : vector giving the fixed parameters. (Must be empty or 
%                 have the same length as the number of parameters. 
%                 Non-fixed parameters must then be given as NaN's)
%                 (default [])
%       .alpha    : Confidence coefficent             (default 0.05)
%       .optimset : optimset structure defining performance of the
%                   optimization routine (see optimset for details))
%
% Example:
%  sz = [1000,1];
%  R1 = rndweib(1,2,0,sz);
%  R2 = rndray(2,sz);
%  opts = fitmarg2d('defaults');
%  opts = parseoptions(opts,'distribution',{'pdfweib','pdfray'});
%  Phat2 = fitmarg2d(R1,R2,opts);
% 
% See also  rndmarg2d,  pdfmarg2d, cdfmarg2d, cdfmargcnd2d

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


% References:
% Plackett, R. L. (1965) "A class of bivariate distributions."
%                                J. Am. Stat. Assoc. 60. 516-22
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 p. 133-134.


%  tested on: matlab 5.2
% history
% Revised pab nov2004
%  -replaced fmins with fminsearch  
% revised pab Dec2003
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vector to structure
% Per A. Brodtkorb 29.01.1999

options = struct('distribution',[],'method','SML','fixpar',[],'start',[],'alpha',0.05,'optimset',optimset('disp','off'));
if nargin==1 && nargout <= 1 && isequal(V,'defaults')
  phato = options;
  return
end
error(nargchk(3,inf,nargin))
options = parseoptions(options,varargin{:});
 
V = V(:);
H = H(:);
[n, m] = size(V);
[n2, m2] = size(H);
if n~=n2
  error('V and H  must have equal size')
end
try
  VDIST=getdistname(lower(options.distribution{1}));
  HDIST=getdistname(lower(options.distribution{2}));
catch
  error('Not enough inputs (distributions missing)')
end


fitvh = str2func(['fit',VDIST]);
fithh = str2func(['fit',HDIST]);

mmethod = options.method;
options.method(1) = '';
ph = fithh(H,options);
pv = fitvh(V,options);

Np1 = numel(pv.params);
Np2 = numel(ph.params);
pdfoptions = struct('numpar',[Np1,Np2],'distribution',[],'condon',0);
pdfoptions.distribution = options.distribution;
if isempty(options.start)
  rho= corrcoef(V,H);
  rho=rho(2,1);
  psi=findPsi(rho);

  phat0=[pv.params ph.params psi ];
else
  phat0 = options.start;
end
if strncmpi(mmethod,'s',1);
  dosearch = true;
else
  dosearch = false;
end
if isempty(options.fixpar)
  if isempty(pv.fixpar)
    pv.fixpar = repmat(nan,size(pv.params));
  end
  if isempty(ph.fixpar)
    ph.fixpar = repmat(nan,size(ph.params));
  end
   options.fixpar = [pv.fixpar,ph.fixpar,nan];
end
if numel(options.fixpar)~=numel(phat0)
  error('Make sure number of variables are consistent')
end
phat0(isfinite(options.fixpar))= [];
phato = mlest(@pdfmarg2d,phat0,{V(:),H(:)},{options,'search',dosearch},pdfoptions);
phato.dataname = inputname(1);
phato.pdfoptions = pdfoptions;
end

function k=findPsi(rho)
  % finds psi by a simple iteration scheme 
  rtol=1e-6; %relativ tolerance
  switch rho
  case 0,  k=1; return
  case 1, k=inf;return
  case -1, k=0; return
  otherwise, 
    if rho<0,
      kold=0.5;
    else
     kold=4;
    end
  end
  kold2=kold-0.001;
  rho2=getrho(kold2);
  for ix=1:500,
    rho1=getrho(kold);
    k=kold-(rho1-rho).*(kold-kold2)./( rho1-rho2);
   %    disp(['k=' num2str(k) ]), pause
    if abs(kold-k)<abs(rtol.*k),break
    elseif k<=0, k=linspace(0,20,1000);
      [tmp I]=min(abs(getrho(k)-rho)); %linear search
      k=k(I);
      break;
    end
    kold=k;
    if ix>=500, disp('could not find k'),end
  end
end

function r=getrho(psi)
r=zeros(size(psi));
k1=find(psi==0);
if any(k1),
 r(k1)=-1;
end
k3=find((psi<1.*psi>0)|psi>1);
if any(k3),
 r(k3)=(psi(k3).^2-1-2.*psi(k3).*log(psi(k3)))./((psi(k3)-1).^2) ;
end
end

function [pvhat,pvci]=distfit(tmp,dist,method,alpha)
  
  if strcmpi(method(2:3),'ml'),
    switch lower(dist(1:2)),
      case 'ra', [pvhat tmp2] = fitray(tmp,0) ;
      case 'tg', [pvhat tmp2] = fitgumb(tmp,0);
      case 'gu', [pvhat tmp2] = fitgumb(tmp,0);
      case 'lo', [pvhat tmp2] = fitnorm(log(tmp),0);
      case 'ga', [pvhat tmp2] = fitgam(tmp,0);
      case 'we', [pvhat tmp2] = fitweib(tmp,0);
      otherwise, error('Unknown distribution')
    end
    if numel(tmp2)~=length(tmp2);
      sa=diag(tmp2)';
    else
      sa = tmp2(:)';
    end
    p_int=[alpha/2; 1-alpha/2];
    pvci = invnorm(repmat(p_int,1,length(pvhat)),[pvhat; pvhat],[sa;sa]);
  else  % MOM fit
    
    error('MOM is not implemented for any distribution')	
%     switch dist(1:2)
%       case {'tg','gu'} ,  
%       case 'we', 
%       case 'lo', 
%       case 'ga', error('MOM is not implemented for Gamma distribution')	 
%       case 'ra', error('MOM is not implemented for Rayleigh distribution')	 
%       otherwise , error('Unknown distribution')
%     end
  end
end

