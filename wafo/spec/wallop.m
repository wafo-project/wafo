function S1=wallop(w1,sdata,plotflag)
%WALLOP  Calculates (and plots) a Wallop spectral density.
% 
% CALL:  S = wallop(w,data,plotflag); 
%        S = wallop(wc,data,plotflag);
%
%        S    = a struct containing the spectral density, see datastructures.
%        w    = angular frequency (default linspace(0,3,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp N]
%               Hm0 = significant wave height (default 7 (m))
%               Tp  = peak period (default 11 (sec))
%               N   = shape factor, i.e. slope for the high frequency
%                     part (default depending on Hm0 and Tp, see below)
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
%  The WALLOP spectrum parameterization used is 
%
%    S(w) = A * G0 * wn^(-N)*exp(-N/(4*wn^4))
%  where 
%    G0 = Normalizing factor related to Bretschneider form
%    A  = (Hm0/4)^2 / wp     (Normalization factor)
%    wn = w/wp       
%    wp = 2*pi/Tp,  angular peak frequency
%    N  = abs((log(2*pi^2)+2*log(Hm0/4)-2*log(Lp))/log(2));
%    Lp = wave length corresponding to the peak frequency, wp.
%
%  If N=5 it becomes the same as the JONSWAP spectrum with 
%  peak enhancement factor gamma=1 or the Bretschneider 
%  (Pierson-Moskowitz) spectrum. 
%
% Example: 
%   S = wallop(1.1,[6.5 10]), plotspec(S)
%  
% See also  jonswap, torsethaugen, simpson

% References:
% Huang, N.E., Long, S.R., Tung, C.C, Yuen, Y. and Bilven, L.F. (1981)
% "A unified two parameter wave spectral model for a generous sea state"
% J. Fluid Mechanics, Vol.112, pp 203-224

% Tested on: matlab 6.0, 5.3
% History:
% revised pab
% - replaced code with a call to mkbretschneider
% revised jr 03.04.2001
%  - added wc to input 
%  - updated information
% revised pab 18.02.2000
%  - normalization so that int S(w) dw = m0
% revised pab 24.01.2000
%  - updated note to 'Wallop Hm0='....
% by pab 01.12.99

monitor=0;

if nargin<3||isempty(plotflag)
  plotflag=0;
end

w = [];
N=[];
if nargin<2||isempty(sdata)
  sdata=[7 11];
else
  switch length(sdata)
    case 1, sdata=[sdata 11];
    case 3, N=sdata(3); sdata=sdata(1:2);
  end
end %

if nargin<1||isempty(w1), 
  wc = 33/sdata(2);
elseif length(w1)==1,   
  wc = w1; 
else
  w = w1;
end
nw = 257;

if isempty(w), w = linspace(0,wc,nw).'; end
 

Hm0=sdata(1);
Tp=sdata(2);
wp=2*pi/Tp;

if monitor
  disp(['Hm0, Tp      = ' num2str([Hm0 Tp])])
end

if isempty(N),
  kp=w2k(wp,0,inf); % wavenumber at peak frequency
  Lp=2*pi/kp; % wave length at the peak frequency
  N=abs((log(2*pi^2)+2*log(Hm0/4)-2*log(Lp))/log(2));
end

% Call modified bretschneider
bspec = mkbretschneider('Hm0',Hm0,'Tp',Tp,'N',N);

S1 = createspec;
S1.S = bspec(w);
S1.w = w;
S1.norm=0; % The spectrum is not normalized


S1.note=['Wallop, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp)];


%Bw=0.06238*M^((M-1)/4)/(4^((M-5)/4)*gamma((M-1)/4))*(1+0.7458*(M+2)^(-1.057));

% Bw = N^((N-1)/4)/(4^((N-5)/4)*gamma((N-1)/4))/16;
% 
% % for w>0 % avoid division by zero
% k=find(w>0);
% S1.S(k)=Bw*Hm0^2/wp*(wp./w(k)).^N.*exp(-N/4*(wp./w(k)).^4);
% 
if plotflag
  plotspec(S1,plotflag)
end

