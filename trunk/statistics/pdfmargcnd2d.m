function y = pdfmargcnd2d(v1,h1,phat)
%PDFMARGCND2D Joint 2D PDF computed as f(x1|X2=x2)*f(x2)
%
% CALL:  f = pdfmargcnd2d(x1,x2,phat) 
%
%      f  = PDF evaluated at x1,x2
%  x1,x2  = evaluation points
%    phat = structure array containing
%           x    = cellarray of distribution parameters
%           dist = cellarray of strings defining the distributions of 
%                  X2 and X1 given X2, respectively. Options are:
%                  'tgumbel', 'gumbel', 'lognormal','rayleigh','weibull',
%                  and 'gamma'.
%
% PDFMARGCND2D evaluates f{x1|X2=x2}*f{x2}. 
%  The parameter(s) of the unconditional distribution of X2,
%  f{x2}, must be in in phat.x{2}. The parameters of the conditional
%  distribution of X1 given X2 must be in phat.x{1}. The first column
%  in phat.x{1} contains the X2 values the parameters in column 2 and 3 are
%  conditioned on.  
% 
% The size of f is the common size of X1 and X2.  
%
% Example: 2D Rayleigh, ie, f(x1)*f(x2)
%   x1=linspace(0,10)';
%   phat0.x={[x1,2*ones(size(x1))] 2 };
%   phat0.dist={'rayl','rayl'};
%   pdfmargcnd2d(2,1,phat0)
%
% See also  fitmargcnd2d rndmargcnd2d pdfmargcnd2d cdfmargcnd2d

%tested on: matlab 5.2
% history:
% revised pab 19.01.2001
% revised pab 03.12.2000
% added truncated weibull and truncated raleigh 
% revised pab 12.11.2000
%  - added ggampdf option
% revised pab 08.02.2000
%   fixed a bug CV -> Cv
%  Per A. Brodtkorb 28.10.98

error(nargchk(3,3,nargin))

dist=phat.dist;

V=v1; 
H=h1; 

y = zeros(max([size(V) ;size(H)]));

if strcmpi('gu', dist{1}(1:2)),
  if strcmpi('gu', dist{2}(1:2)),
    k=find(H>-inf);
  else
    k=find(H>0);
  end
elseif strcmpi('gu', dist{2}(1:2)),
  k = find(V > 0 );
else
  k = find(V > 0 & H>0);
end

PH=phat.x{2};


if any(k),     
    switch lower(dist{2}(1:2))
      case 'tr' ,  pdf1=pdfraymod(H(k),PH(1),PH(2));
      case 'ra' ,  pdf1=pdfray(H(k),PH);
      case 'we' ,  pdf1=pdfweib(H(k),PH(1),PH(2));
      case 'tw' ,  pdf1=pdfweibmod(H(k),PH(1),PH(2),PH(3));
      case 'gu' ,  pdf1=pdfgumb(H(k),PH(1),PH(2),0);
      case 'tg' ,  pdf1=pdfgumb(H(k),PH(1),PH(2),1);
      case 'ga' ,  pdf1=pdfgam(H(k),PH(1),PH(2));
      case 'gg' ,  pdf1=pdfgengam(H(k),PH(1),PH(2),PH(3));
      case 'lo' ,  pdf1=pdflognorm(H(k),PH(1),PH(2));
      otherwise, error('unknown distribution')
    end 
    [Av , Bv, Cv]=margcnd2dsmfun(phat,H(k)); %parameters of V given H 
   switch lower(dist{1}(1:2))
     case 'tr', y(k) =  pdf1.*pdfraymod(V(k),Av,Bv);
     case 'ra', y(k) =  pdf1.*pdfray(V(k)-Cv,Av);
     case 'gu', y(k) =  pdf1.*pdfgumb(V(k)-Cv,Av,Bv,0);
     case 'tg', y(k) =  pdf1.*pdfgumb(V(k)-Cv,Av,Bv,1);
     case 'lo', y(k) =  pdf1.*pdflognorm(V(k)-Cv,Av,Bv);
     case 'ga', y(k) =  pdf1.*pdfgam(V(k)-Cv,Av,Bv);	
     case 'gg', y(k) =  pdf1.*pdfgengam(V(k),Av,Bv,Cv);	
     case 'we', y(k) =  pdf1.*pdfweib(V(k)-Cv,Av,Bv);
     case 'tw', y(k) =  pdf1.*pdfweibmod(V(k),Av,Bv,Cv);
     otherwise, error('Unknown distribution')
   end
end

%y(find(isnan(y)|isinf(y)))=0;

