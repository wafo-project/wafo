function p = bmr00pdf2(Hd,Scf,Hs,Tz,normalizedInput)
%BMR00PDF2 Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
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
%    The size of f is the common size of  the input arguments
%   
% Example:
%  Hs = 7;Tz=10;  
%  h = linspace(0,3*Hs)'; 
%  s = linspace(0,3*5/4*Hs/Tz^2)';
%  f = bmr00pdf2(h,s,Hs,Tz);
%  pdfplot(f)
%
% See also mk87pdf2

% Reference 
%  Brodtkorb, P.A. and Myrhaug, D. and Rue, H. (2000)
%  "Joint Distributions of Wave Height and Wave Steepness Parameters",
% In Proc. 27'th Int. Conf. on Coastal Eng., ICCE, Sydney, Australia },
%  vol. 1, pp. 545--558, Paper No. 162

% By pab 15 July 2004

%error(nargchk(3,5,nargin))
narginchk(3,5)
if nargin < 5||isempty(normalizedInput), normalizedInput = 0; end
if nargin < 4||isempty(Tz),  Tz = 8; end
if nargin < 3||isempty(Hs),  Hs = 6; end

p     = createpdf(2);

[X,Y] = meshgrid(Scf,Hd);
p.f   = bmr00pdf(Y,X,Hs,Tz,0,normalizedInput);
p.x{1}=Scf(:);
p.x{2}=Hd(:);
p.labx{1}='Crest front steepness';
p.labx{2}='Wave height [m]';

p.title='Brodtkorb et.al. (2000) density of crest front steepness and wave height.';
  
[p.cl p.pl]=qlevels(p.f,[10:20:90 95 99 99.9]);
