function S = oscspec(sdata,z,b,s)
% OSCSPEC Spectral density for a harmonic oscillator
%         driven by white noise
%  
% CALL:  S = oscspec(sdata,z,b,s);
%
%        S     = the spectral density (structure array)
% 	 sdata = the data vector [wl wu n], where 
% 
%	    wl = lower truncation frequency  (default 4/257)
%	    wu = upper truncation frequency  (default 4)
%	    n  = number of evaluation points (default 257)
%
%       z,b,s = parameters in the equation (eq. 1) for the oscillator.
%               (default z=0.1, b=1, s=1)
%        
% Let W  be white noise, then the oscillator X is defined by
%  
%     X''(t) + 2bz X'(t) + b^2 X(t) = s W(t)    
%
% The angular peak frequency is given by wp = b/sqrt(1-4*z^2).
% Important parameter values:
%    0<z<1, b=1, s = 2*sqrt(z) : Normalized linear oscillator,  Var(X(t))=Var(X'(t))=1.
%   
% Example: 
%  data = [0.01 4 275];
%  S    = oscspec(data,[],2.5);  % Peak frequency at w=2.5
%
% See also duffsim

% Tested on: Matlab 5.3
% History: 
% -Revised pab Feb2007
% - changed numenclature: w->b in order to avoid confusion
% - multiplied spectrum with 2 to get correct results.
% Correction by PJ 07-Jul-2005
%   Changed 'break' to 'return'
% revised es 25.05 00 small modifications of help text   
% Modified by jr 14.01.2000
% - updated check of nargins
% Modified by jr 12.01.2000
% - new names of variables and parameters
% - check of nargins introduced
% - updated documentation
% Modified by ir 11.01.2000
% - structure array introduced
% - parameters w,s allowed as input
% By Mats Frendahl 1993
  
if nargin<1 || isempty(sdata), sdata=[4/257 4 257];end
if nargin<2 || isempty(z), z=0.1;end
if nargin<3 || isempty(b), b=1;end
if nargin<4 || isempty(s), s=1;end

if (z<0)||(z>1)
  error('WAFO:OSCSPEC',' The parameter  z  must be in [0,1]. Program will terminate.')
  %return
end

wl = sdata(1); 
wu = sdata(2); 
n  = sdata(3); 

w=linspace(0,wu,n);
if b~=0
  wn2 = (w/b).^2; 
  spv=(s/b)^2/(pi)./((1-wn2).^2+(2*z)^2*wn2);
else
  spv = zeros(size(w));
  spv(w>0)=s^2./(pi*w(w>0).^4);
end

S=createspec;
S.S=spv; 
S.w=w;
S.type='freq';
S.note=sprintf('Spectrum, harmonic oscillator, wp = %2.1f',b/sqrt(1-4*z^2));
S.S(w<wl)=0;
S.S(1)=0; % must be zero at zero freq since discrete spectrum
%S=floor(S*1e5+.5)/1e5;



return

