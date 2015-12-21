function p = mk87pdf2(Hd,Scf,Hs,Tz)
%MK87PDF2 Myrhaug and Kjeldsen (1987) joint (Scf,Hd) PDF. 
%
% CALL:  f = mk87pdf2(Hd,Scf,Hs,Tz)
%
%    f  = pdf structure evaluated at meshgrid(Scf,Hd) 
%   Hd  = zero down crossing wave height
%   Scf = crest front steepness
%   Hs  = significant wave height.
%   Tz  = average zero down crossing period.  
%
% MK87PDF2 returns the joint PDF of (Scf, Hd) given Hs and Tz,
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, given
% the seastate. The root mean square values of Hd and Scf (Hrms,Erms) are
% related to the significant waveheight and the average zero down
% crossing period by:
%             Hrms = 0.715*Hs;
%             Erms = 0.0202+0.826*Hs/(Tz^2);
%
% Example:
%  Hs = 7;Tz=10;  
%  h = linspace(0,3*Hs)'; 
%  s = linspace(0,3*Hs/Tz^2)';
%  f = mk87pdf2(h,s,Hs,Tz);
%  pdfplot(f)
%
% See also  mk87pdf, createpdf

%   References:
%   Myrhaug, D. and Kjelsen S.P. (1987) 
%  'Prediction of occurences of steep and high waves in deep water'.
%   Journal of waterway, Port, Coastal and Ocean Engineers, Vol. 113, pp 122--138
%
%   Myrhaug & Dahle (1984) 
%   'Parametric modelling of joint probability density 
%   distributions for steepness and asymmetry in deep water waves'
%

% tested on: matlab 5.1
% history:
% revised pab 09.08.2003
% Changed input + updated help header  
% revised pab 01.04.2001
% -added example
% revised pab 08.02.2000
%  - 
% by  Per A. Brodtkorb 1998

error(nargchk(3,4,nargin))

if nargin < 4||isempty(Tz),  Tz = 8; end
if nargin < 3||isempty(Hs),  Hs = 6; end

p     = createpdf(2);

[X,Y] = meshgrid(Scf,Hd);
p.f   = mk87pdf(Y,X,Hs,Tz);
p.x{1}=Scf(:);
p.x{2}=Hd(:);
p.labx{1}='Crest front steepness';
p.labx{2}='Wave height [m]';

p.title='Myrhaug and Kjeldsen (1987) density of crest front steepness and wave height.';
  
[p.cl p.pl]=qlevels(p.f,[10:20:90 95 99 99.9]);
