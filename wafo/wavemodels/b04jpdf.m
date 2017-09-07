function f = b04jpdf(Hd,Scf,Hm0,Tm02, condon,normalizedInput)
%B04JPDF Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
%   
%    CALL:  f = b04jpdf(Hd,Scf,Hm0,Tm02)
%   
%       f  = density 
%      Hd  = zero down crossing wave height
%      Scf = crest front steepness
%      Hm0 = significant wave height.
%     Tm02 = average zero down crossing period.  
%   
%    B04JPDF returns the joint PDF of (Scf, Hd) given Hm0 and Tm02,
%    i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%    given the seastate. The root mean square values of Hd and Scf 
%    (Hrms,Erms) are related to the significant waveheight and the
%    average zero down crossing period by:
%                Hrms = Hm0/sqrt(2);
%                Erms = 5/4*Hm0/(Tm02^2);
%    This distribution  is fitted to storm waves from the Japan Sea.
%    The size of f is the common size of  the input arguments
%   
%    Example:
%     Hs = 7;Tz=10;  
%     h  = linspace(0,3*Hs)'; 
%     s  = linspace(0,4*Hs/Tz^2)';
%     [S,H] = meshgrid(s,h);
%     contour(s,h,b04jpdf(H,S,Hs,Tz))
%
% See also mk87pdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.

% By pab 15 July 2004
 

 
 
%error(nargchk(3,6,nargin))
narginchk(3,6)  
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
  % Error 0.07 for cAx and 0.01 for cBx
  cA1 = [   2.42767540920461,...
	    -5.12101444148408,...
	    4.49060259316430];
  
  cA2 =[   0.00572707656955,...
	   -0.10679331978077,...
	   0.81081632991506,...
	   -3.19389926947809,...
	   6.81410759646213,...
	   -7.14860241955959,...
	   2.40890619464040,...
	   1.00000000000000];
  cB1 =[  -0.70310486398236,...
	  1.74679333441188,...
	  -2.47501962948775];
  
  cB2 =[   0.00536662963627,...
	   -0.10086105895100,...
	   0.77434723641863,...
	   -3.10647077544176,...
	   6.88647211568827,...
	   -7.71847113553527,...
	   3.19755172037110,...
	   1.00000000000000];

end
A1 = polyval(cA1,h)./polyval(cA2,h);
B1 = exp(polyval(cB1,h)./polyval(cB2,h));
h0 = 4.49999999999989;
k = find(h>h0);
if any(k)
  % Linear extrapolation
 slope1 = 1.53878793498107;
 
 slope2 =  0.01549012000025;
     
  A1(k) =  (h(k)-h0)*slope1 + 17.80538369601295;
  B1(k) = (h(k)-h0)*slope2  + 0.15016765927218;
end
% Rayleigh parameter
B0 = 0.67768680795379; 


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
% $$$ m = 5;
% $$$ k = m+1;
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
% $$$ 
% $$$ 
