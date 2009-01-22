function [F ,err] = cdfweib2d(x1,x2,varargin)
%CDFWEIB2D Joint 2D Weibull cumulative distribution function  
%
%  CALL: [F,err] = cdfweib2d(x1,x2,A1,B1,A2,B2,C12,options)
%        F = cdfweib2d(x1,x2,phat,options)
%     F    = joint CDF, i.e., Prob(X1<x1 X2<x2)
%    err   = estimated error
%    x1,x2 = coordinates
%   A1, A2 = scale parameters    
%   B1, B2 = shape parameters
%      C12 = interaction parameter between X1 and X2
%     phat = Distribution parameter struct
%            as returned from FITWEIB2D.  
%  options = struct with fieldnames:
%     .logp    : if TRUE, probability, p, returned as log(p).
%     .condon  : 0 it returns the the regular cdf of X1 and X2 (default)
%                1 it returns the conditional cdf given X1
%                2 it returns the conditional cdf given X2
%     .releps  : specified relative tolerance (Default 1e-3) 
%
% The size of F is the common size of X1 and X2.
% The CDF is defined by its PDF:
%
%  f(X1,X2)=B1*B2*xn1^(B1-1)*xn2^(B2-1)/A1/B1/N*...
%           exp{-[xn1^B1 +xn2^B2 ]/N }*I0(2*C12*xn1^(B1/2)*xn2^(B2/2)/N) 
%   where 
%      N=1-C12^2, xn1=X1/A1,  xn2=X2/A2 and 
%      I0 is the modified  bessel function of zeroth order.
%
%  (The marginal distribution of X1 is weibull with parameters A1 and B1) 
%  (The marginal distribution of X2 is weibull with parameters A2 and B2) 
%
% Example:
%   x = linspace(0,6,200); [X1,X2] = meshgrid(x); 
%   phat = {2 2  3 2.5 .8};
%   f = pdfweib2d(X1,X2,phat{:});
%   contour(x,x,f), hold on,
%   plot( [0  2  2 0 0], [0 0 1 1 0],'g-'), hold off % Mark the region
%   cdfweib2d(2,1,phat{:})  % Calculate the probability of marked region
% 
% See also  pdfweib2d, prbweib2d, rndweib2d, fitweib2d, momweib2d, quad2dg, gaussq

% tested on: matlab5.1
%history:
% revised pab Dec2003
%  call weibstat -> momweib
% revised pab 29.10.2000
% - updated to wstats
% - added example
%  Per A. Brodtkorb 17.11.98

% NB! weibpdf must be modified to correspond to
% pdf=x^(b-1)/a^b*exp(-(x/a)^b) 

error(nargchk(3,15,nargin))
Np = 5;
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'condon',0,'releps',1e-3,...
  'meshgrid',false,'wdata',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a1,b1,a2 b2,c12] = deal(params{:});

if isempty(a1)||isempty(b1)||isempty(a2)||isempty(b2)||isempty(c12)
  error('Requires either 7 input arguments or that input argument 3 is FDATA.'); 
end

[icode x1 x2 ] = iscomnsize(x1,x2);
if  ~icode
  error('Requires non-scalar arguments to match in size.');
end
F = zeros(size(x1));
err=F;

if 0
  % This is a trick to get the html documentation correct.
  k = pdfweib2d(1,1,2,3);
end

tol = options.releps;

switch options.condon
  case 0,% regular cdf   
    k = find( x1 > 0 & x2>0);
    if 1, 
       if any(k), 
         [F(k) err(k)]= gaussq2d(@pdfweib2d,0,x1(k),0 , x2(k),  tol,params{:});
       end
     else% divide the integral in to several parts this is not correct yet
%       [m1]=momweib(a1,b1); 
%       [m2]=momweib(a2,b2);
%        x1s=min(m1,x1);
%        x2s=min(m2,x2);
%        if any(k)
%          [F(k) err(k)]=gaussq2d(@pdfweib2d,0,x1s(k),0 , x2s(k),tol/2,params{:},condon);
%        end
%        k1=find( x1(k)>m1& x2(k)<=m2);
%        if any(k1),
%          [tmp1 tmp2]=gaussq2d(@pdfweib2d,x1s(k(k1)),  x1(k(k1)),0,   x2s(k(k1)), tol/4,params{:});
%          F(k(k1)) =F(k(k1))+tmp1;
%          err(k(k1))=err(k(k1))+tmp2;
%        end
%        k1=find( x1(k)<=m1& x2(k)>m2);
%        if any(k1),
%          [tmp1 tmp2]=gaussq2d(@pdfweib2d,0,x1s(k(k1)), x2s(k(k1)),   x2(k(k1)),   tol/4,params{:});
%          F(k(k1)) =F(k(k1)) +tmp1;
%          err(k(k1))=err(k(k1))+tmp2;
%        end
%        k1=find(m1<x1(k)& m2<x2(k));
%        if any(k1),
%          [tmp1 tmp2]=gaussq2d(@pdfweib2d,x1s(k(k1)),x1(k(k1)), x2s(k(k1)),   x2(k(k1)),   tol/4,param,condon);
%          F(k(k)) =F(k(k)) +tmp1;
%          err(k(k))=err(k(k))+tmp2;
%        end
     end
   case 1,%conditional CDF given x1  
     k = find( (x1 > 0) & (x2>0) & (c12(ones(size(x2))) ~=0 ));
     if any(k),
       [F(k) err(k)]=gaussq(@pdfweib2d,0,x2(k),tol,[],x1(k),a2,b2,a1,b1,c12,'condon',2);%param([3:4 1:2 5]),2);   
     end
     k = find( (x1==0)| (c12(ones(size(x2))) ==0));
     if any(k),
       F(k) =cdfweib(x2(k),params{3:4}); 
     end
     %for ix=1:length(x1(:)),      
     %[F(ix) err(ix)]=quadg('pdfweib2d',0,x2(ix),tol,[],x1(ix),param([3:4 1:2 5]),2);      
     %end
   case 2,%conditional CDF given x2  
     k = find( (x1 > 0) & (x2>0) & (c12(ones(size(x2))) ~=0));
     if any(k),
       [F(k) err(k)]=gaussq(@pdfweib2d,0,x1(k),tol,[],x2(k),a1,b1,a2,b2,c12,'condon',2);%param,2); 
     end
     k = find( (x2==0)| (c12(ones(size(x2))) ==0));
      if any(k),
       F(k) =cdfweib(x1(k),params{1:2}); 
     end
     %for ix=1:length(x2(:)), 
      % [F(ix) err(ix)]=quadg('pdfweib2d',0,x1(ix),tol,[],x2(ix),param,2);      
     %end
 end    
 
% make sure 0 <= y<=1 
k2=find(F<0);
if any(k2)
  F(k2)=zeros(size(k2));
end
k2=find(F>1);
if any(k2)
  F(k2)=ones(size(k2));
end
if any(isnan(F)),
  disp(['Warning: there are : ', num2str(sum(isnan(y))),' NaNs']);
end
