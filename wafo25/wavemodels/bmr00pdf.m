function f = bmr00pdf(Hd,Scf,Hm0,Tm02, condon,normalizedInput)
%BMR00PDF Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
%   
%    CALL:  f = bmr00pdf(Hd,Scf,Hm0,Tm02)
%   
%       f  = density 
%      Hd  = zero down crossing wave height
%      Scf = crest front steepness
%      Hm0 = significant wave height.
%     Tm02 = average zero down crossing period.  
%   
%    BMR00PDF returns the joint PDF of (Scf, Hd) given Hm0 and Tm02,
%    i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%    given the seastate. The root mean square values of Hd and Scf 
%    (Hrms,Erms) are related to the significant waveheight and the
%    average zero down crossing period by:
%                Hrms = Hm0/sqrt(2);
%                Erms = 5/4*Hm0/(Tm02^2);
%    This is a revised distribution of MK87 and is fitted to storm waves
%    from 1995 obtained from the Draupner field in the North Sea. 
%      The size of f is the common size of  the input arguments
%   
%    Example:
%     Hs = 7;Tz=10;  
%     h  = linspace(0,3*Hs)'; 
%     s  = linspace(0,4*Hs/Tz^2)';
%     [S,H] = meshgrid(s,h);
%     contour(s,h,bmr00pdf(H,S,Hs,Tz))
%
% See also mk87pdf

% Reference 
%  Brodtkorb, P.A. and Myrhaug, D. and Rue, H. (2000)
%  "Joint Distributions of Wave Height and Wave Steepness Parameters",
% In Proc. 27'th Int. Conf. on Coastal Eng., ICCE, Sydney, Australia },
%  vol. 1, pp. 545--558, Paper No. 162

% By pab 15 July 2004
 

 
 
error(nargchk(3,6,nargin))
  
if nargin < 6||isempty(normalizedInput),  normalizedInput=0; end
if nargin < 5||isempty(condon),  condon=0; end % regular pdf is returned
if nargin < 4||isempty(Tm02),  Tm02=8;end
if nargin < 3||isempty(Hm0),  Hm0=6;end

[csize, Scf,Hd,Tm02,Hm0] = comnsize(Scf,Hd,Tm02,Hm0);
if any(isnan(csize))
    error('Requires non-scalar arguments to match in size.');
end
if (normalizedInput>0)
  Hrms = 1;
  Erms = 1;
else
  %Hm00 = 6.7;
  %Tm020 = 8.3
  Hrms = Hm0/sqrt(2);
  Erms = 5/4*Hm0./(Tm02.^2); 
  %Erms = (0.0202 + 0.826*Hm0./(Tm02.^2))/2; 
end
s = Scf./Erms;
h = Hd./Hrms;

% Conditinal gamma dist.  parameters for steepness
if 1,
  % Error 0.02 for cAx and 0.02 for cBx
  cA1 = [ 4.19193631070257,...
	  -12.03484717258579,...
	  15.95695603911886,...
	  -7.31618437661838,...
	  3.03590705490789];
  cA2 = [  -0.02838644340635,...
	   0.29409707273437,...
	   -0.79138075176131,...
	   1.20024905431293,...
	   -0.58192823600509,...
	   1.00000000000000];
  cB1 =[ -5.96970985003487,...
  23.24982082461128,...
 -37.69215199928053,...
   5.18520245368619,...
  -1.87015799594426];
  cB2 = [ -0.47236526233572,...
	  5.71130418016660,...
	  -17.45884089250511,...
	  21.60781934904390,...
	  0.11061694237412,...
	  1.00000000000000];
else
  % Error 0.08 for cAx and 0.03 for cBx
  cA1 = [3.39863744229230, 0.42853512533872,2.97401909675030];
  cA2 = [0.26749624884949, -1.47019799044842,2.17197357932780,1];
  cB1 = [-15.50443060097295, -0.21474910268201,  -1.87709428385351];
  cB2 = [-0.49720581030185, 6.35119691685354, 3.88968350685379, 1];
 end
A1 = polyval(cA1,h)./polyval(cA2,h);
B1 = exp(polyval(cB1,h)./polyval(cB2,h));
h0 = 3.79937912772197;
k = find(h>h0);
if any(k)
  % Linear extrapolation
 slope1 = 14.03628528937437;
 slope2 = -0.37497306844631;
  A1(k) =  (h(k)-h0)*slope1 + 36.15168757587687;
  B1(k) = exp((h(k)-h0)*slope2 + log(0.05673842727364));
  %3.37291820522257  30.10141068490231   0.06979731246760
  %3.41168737999524  30.65143585680909   0.06861014108633
end
% Rayleigh parameter
B0 = 0.69233682309018;



%f   = zeros(size(h));

switch condon,
 case 0, % regular pdf is returned 
  f = pdfray(h,B0).*pdfgam(s,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfgam(s,A1,B1);
 case 3, % secret option  used by mk87stat: returns x2*p(x2|x1) 
  f = s.*pdfgam(s,A1,B1);
 case 4, % secret option  used by mk87stat: returns x2.^2*p(x2|x1) 
  f = s.^2.*pdfgam(s,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by bmr00cdf
  f = pdfray(h,B0).*cdfgam(s,A1,B1);
 case 6, % P(V|h) is returned special case used by bmr00cdf
  f = cdfgam(s,A1,B1);
 case 7, % p(h)*(1-P(V|h)) is returned special case used by bmr00cdf
  f = pdfray(h,B0).*(1-cdfgam(s,A1,B1));
  otherwise, error('unknown option')
end
if condon~=6
  f = f./Hrms./Erms;
end
return
 
 
 
 % Code used for finding the rational polynomials
% $$$ NFAC=8;
% $$$ BIG=1e30;
% $$$ MAXIT=5;
% $$$ 
% $$$ dev=BIG;
% $$$ 
% $$$ 
% $$$ 
% $$$ m = 7;
% $$$ k = 7;
% $$$ ncof=m+k+1;
% $$$ npt=NFAC*ncof; 
% $$$ % Number of points where function is avaluated, i.e. fineness of mesh
% $$$ 
% $$$ a = 0; b = 2.5
% $$$ 
% $$$ x = zeros(npt,1);
% $$$ ix1=1:floor(npt/2-1);
% $$$ x(ix1)=a+(b-a).*sin(pi/2*(ix1-1)./(npt-1)).^2;
% $$$ ix2=floor(npt/2):npt;
% $$$ x(ix2)=b-(b-a).*sin(pi/2*(npt-ix2)./(npt-1)).^2;
% $$$ [Av , Bv, Cv]=dist2dsmfun(sphat,x);
% $$$ 
% $$$ 
% $$$ [cA1, cA2] = ratlsq([x,Av],m,k,a,b);
% $$$ 
% $$$ [cB1, cB2] = ratlsq([x,log(Bv)],m,k,a,b);
% $$$ %[cA1, cA2] = ratlsq([x,Av],a,b,m,k);
% $$$ subplot(2,1,1)
% $$$  plot(x,polyval(cA1,x)./polyval(cA2,x),'g',x,Av,'r')
% $$$  subplot(2,1,2)
% $$$  plot(x,exp(polyval(cB1,x)./polyval(cB2,x)),'g',x,Bv,'r')
% $$$  
% $$$  
% $$$  
% $$$ m = 3;
% $$$ k = m;
% $$$ ncof=m+k+1;
% $$$ npt=NFAC*ncof; 
% $$$ % Number of points where function is avaluated, i.e. fineness of mesh
% $$$ 
% $$$ a = 0; b = 2
% $$$ 
% $$$ x = zeros(npt,1);
% $$$ ix1=1:floor(npt/2-1);
% $$$ x(ix1)=a+(b-a).*sin(pi/2*(ix1-1)./(npt-1)).^2;
% $$$ ix2=floor(npt/2):npt;
% $$$ x(ix2)=b-(b-a).*sin(pi/2*(npt-ix2)./(npt-1)).^2;
% $$$ [Av , Bv, Cv]=dist2dsmfun(sphat,x);
% $$$ 
% $$$ 
% $$$ [cA1, cA2] = ratlsq([x,Av],m,k,a,b);
% $$$ 
% $$$ [cB1, cB2] = ratlsq([x,log(Bv)],m,k,a,b);
% $$$ %[cA1, cA2] = ratlsq([x,Av],a,b,m,k);
% $$$ subplot(2,1,1)
% $$$  plot(x,polyval(cA1,x)./polyval(cA2,x),'g',x,Av,'r')
% $$$  subplot(2,1,2)
% $$$  plot(x,exp(polyval(cB1,x)./polyval(cB2,x)),'g',x,Bv,'r')
% $$$  
% $$$  N = 12
% $$$  x1 = chebroot(N).'*(b-a)/2+(b+a)/2 ;
% $$$  [Av , Bv, Cv]=dist2dsmfun(sphat,x1);
% $$$  cA1 =chebfit([x1 Av],n,a,b);
% $$$  cB1 =chebfit([x1 Bv],n,a,b);
% $$$  
% $$$  dist2dparamplot(phat,sphat)
% $$$  subplot(2,1,1), hold on
% $$$  plot(x,chebval(x,cA1,a,b),'g')
% $$$  subplot(2,1,2), hold on
% $$$  plot(x,chebval(x,cB1,a,b),'g')