function p = b04jpdf2(Hd,Scf,Hs,Tz,normalizedInput)
%B04JPDF2 Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
%   
%    CALL:  f = b04jpdf2(Hd,Scf,Hm0,Tm02)
%   
%       f  = density 
%      Hd  = zero down crossing wave height
%      Scf = crest front steepness
%      Hm0 = significant wave height.
%     Tm02 = average zero down crossing period.  
%   
%    B04JPDF2 returns the joint PDF of (Scf, Hd) given Hm0 and Tm02,
%    i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%    given the seastate. The root mean square values of Hd and Scf 
%    (Hrms,Erms) are related to the significant waveheight and the
%    average zero down crossing period by:
%                Hrms = Hm0/sqrt(2);
%                Erms = 5/4*Hm0/(Tm02^2);
%    This distribution  is fitted to storm waves from Japan Sea.
%    The size of f is the common size of  the input arguments
%   
% Example:
%  Hs = 7;Tz=10;  
%  h = linspace(0,3*Hs)'; 
%  s = linspace(0,3*5/4*Hs/Tz^2)';
%  f = b04jpdf2(h,s,Hs,Tz);
%  pdfplot(f)
%
% See also mk87pdf2

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
  

% By pab 15 July 2004

error(nargchk(3,5,nargin))
if nargin < 5||isempty(normalizedInput), normalizedInput = 0; end
if nargin < 4||isempty(Tz),  Tz = 8; end
if nargin < 3||isempty(Hs),  Hs = 6; end

p     = createpdf(2);

[X,Y] = meshgrid(Scf,Hd);
p.f   = b04jpdf(Y,X,Hs,Tz,0,normalizedInput);
p.x{1} = Scf(:);
p.x{2} = Hd(:);
%p.labx{1}='Crest front steepness';
%p.labx{2}='Wave height [m]';
if (normalizedInput)
  p.labx={'Scf', 'Hd'};
  p.norm = 1;
else
  p.norm=0;
  p.labx={'Scf', 'Hd [m]'};
end

p.title='Brodtkorb (2004) density of (Scf, Hd) Japan Sea.';
  
[p.cl p.pl]=qlevels(p.f,[10:20:90 95 99 99.9]);
