function p = mk87pdf(Hd,Scf,Hs,Tz, condon,normalizedInput)
%MK87PDF Myrhaug and Kjeldsen (1987) joint (Scf,Hd) PDF. 
%
% CALL:  f = mk87pdf(Hd,Scf,Hs,Tz)
%
%    f  = density 
%   Hd  = zero down crossing wave height
%   Scf = crest front steepness
%   Hs  = significant wave height.
%   Tz  = average zero down crossing period.  
%
% MK87PDF returns the joint PDF of (Scf, Hd) given Hs and Tz,
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, given
% the seastate. The root mean square values of Hd and Scf (Hrms,Erms) are
% related to the significant waveheight and the average zero down
% crossing period by:
%             Hrms = 0.715*Hs;
%             Erms = 0.0202+0.826*Hs/(Tz^2);
%
%   The size of f is the common size of  the input arguments
%
% Example:
%  Hs = 7;Tz=10;  
%  h  = linspace(0,3*Hs)'; 
%  s  = linspace(0,4*Hs/Tz^2)';
%  [S,H] = meshgrid(s,h);
%  contour(s,h,mk87pdf(H,S,Hs,Tz))
%
% See also  mk87pdf2, pdfweib, pdflognorm

%   References:
%   Myrhaug, D. and Kjelsen S.P. (1987) 
%  'Prediction of occurences of steep and high waves in deep water'.
%   Journal of waterway, Port, Coastal and Ocean Engineers, 
%   Vol. 113, pp 122--138
  
%
%   Myrhaug & Dahle (1984) 
%   'Parametric modelling of joint probability density 
%   distributions for steepness and asymmetry in deep water waves'
%

% tested on: matlab 5.1
% history:
% revised pab 09.08.2003  
% Changed input + updated help header
% revised pab 23.01.2001
% - no longer dependent on stats toolbox only wstats
% revised pab 10.02.2000
% - fixed normalization
% by  Per A. Brodtkorb 1998


error(nargchk(3,6,nargin))
  
if nargin < 6||isempty(normalizedInput),  normalizedInput=0; end
if nargin < 5||isempty(condon),  condon=0; end % regular pdf is returned
if nargin < 4||isempty(Tz),  Tz=8;end
if nargin < 3||isempty(Hs),  Hs=6;end

[icode, Scf,Hd,Tz,Hs] = iscomnsize(Scf,Hd,Tz,Hs);
if ~icode 
    error('Requires non-scalar arguments to match in size.');
end
if (normalizedInput>0)
  Hrms = 1;
  Erms = 1;
else
  Hrms = 0.715*Hs;
  Erms = 0.0202 + 0.826*Hs./(Tz.^2); 
end

s = Scf./Erms;
h = Hd./Hrms;

sig = (-0.21*atan(2*(h-1.4))+0.325);
sig = max(sig,eps); % pab fix: make sure it is positive
%p   = zeros(size(h));

% NB! weibpdf must be modified to correspond to
% pdf=x^(b-1)/a^b*exp(-(x/a)^b) or else insert
% weibpdf=2.39.*h.^1.39/(1.05^2.39).*exp(-(h./1.05).^2.39);
switch condon,
 case 0, % regular pdf is returned 
  p = pdfweib(h,1.05,2.39).*pdflognorm(s,my(h),sig);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  p = pdflognorm(s,my(h),sig);
 case 3, % secret option  used by mk87stat: returns x2*p(x2|x1) 
  p = s.*pdflognorm(s,my(h),sig);
 case 4, % secret option  used by mk87stat: returns x2.^2*p(x2|x1) 
  p = s.^2.*pdflognorm(s,my(h),sig);
 case 5, % p(h)*P(V|h) is returned special case used by mk87cdf
  p = pdfweib(h,1.05,2.39).*cdflognorm(s,my(h),sig);
 case 6, % P(V|h) is returned special case used by mk87cdf
  p = cdflognorm(s,my(h),sig);
 case 7, % p(h)*(1-P(V|h)) is returned special case used by mk87cdf
  p = pdfweib(h,1.05,2.39).*(1-cdflognorm(s,my(h),sig));
  otherwise, error('unknown option')
end
if condon~=6
  p = p./Hrms./Erms;
end
return

function y = my(h)
y   = zeros(size(h));
ind = (h <= 1.7);
h1  = h(ind);
y(ind) = 0.024-1.065.*h1+0.585.*h1.^2;

%ind=find(h > 1.7);
%h2=h(~ind);
y(~ind) = 0.32*atan(3.14*(h(~ind)-1.7)) - 0.096;
return
