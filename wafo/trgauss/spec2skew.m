function [skew, kurt, ma, sa, Hs ,Hd]=spec2skew(S,h,method)
%SPEC2SKEW Estimates the moments of 2'nd order non-linear waves 
%
% CALL: [skew, kurt, mean, sigma] = spec2skew(S,h,method);
%
%   skew, kurt, 
%   mean, sigma = skewness, kurtosis, mean and standard deviation,
%                 respectively, of 2'nd order waves to the leading
%                 order. (skew=kurt=0 for a Gaussian process)
%             S = spectral density structure
%             h = water depth (default S.h) 
%   	 method = 'approximate' method due to Marthinsen & Winterstein (default)
%                 'eigenvalue'  method due to Kac and Siegert
%
%  Skewness = kurtosis-3 = 0 for a Gaussian process.
%  The mean, sigma, skewness and kurtosis are determined as follows:
%
%  method == 'approximate':  due to Marthinsen and Winterstein
%    mean  = 2 * int Hd(w1,w1)*S(w1) dw1
%    sigma = sqrt(int S(w1) dw1)
%    skew  = 6 * int int [Hs(w1,w2)+Hd(w1,w2)]*S(w1)*S(w2) dw1*dw2/m0^(3/2)
%    kurt  = (4*skew/3)^2
%
%  where Hs = sum frequency effects  and Hd = difference frequency effects
%
% method == 'eigenvalue'
%  
%   mean  = sum(E);                    
%   sigma = sqrt(sum(C^2)+2*sum(E^2));   
%   skew  = sum((6*C^2+8*E^2).*E)/sigma^3;
%   kurt  = 3+48*sum((C^2+E^2).*E^2)/sigma^4;
%
%  where
%   h1 = sqrt(S*dw/2);
%   C  = (ctranspose(V)*[h1;h1]);
%   and E and V are the eigenvalues and eigenvectors, respectively, of the 2'order 
%   transfer matrix. S is the spectrum and dw is the frequency spacing of S.
%
% Example:  Simulate a Transformed Gaussian process:
%  Hm0=7;Tp=11;
%  S = jonswap([],[Hm0 Tp]); [sk, ku, me]=spec2skew(S);
%  g=hermitetr([],[Hm0/4 sk ku me]);  g2=[g(:,1), g(:,2)*Hm0/4];
%  ys = spec2sdat(S,15000);   % Simulated in the Gaussian world
%  xs = gaus2dat(ys,g2);      % Transformed to the real world
%
% See also  hermitetr, ochitr, lc2tr, dat2tr

% References:
% Langley, RS (1987)
% 'A statistical analysis of nonlinear random waves'
% Ocean Engineering, Vol 14, No 5, pp 389-407
%
% Marthinsen, T. and Winterstein, S.R (1992)
% 'On the skewness of random surface waves'
% In proceedings of the 2nd ISOPE Conference, San Francisco, 14-19 june.
%
% Winterstein, S.R, Ude, T.C. and Kleiven, G. (1994)
% "Springing and slow drift responses:
%  predicted extremes and fatigue vs. simulation"
% In Proc. 7th International behaviour of Offshore structures, (BOSS)
% Vol. 3, pp.1-15

% tested on: matlab 5.2 
% History
% revised pab 22.07.2002
%  -fixed a bug in the calculaton of the mean for the 2'nd order process.
%  - added method
%  -updated the help header accordingly
% revised pab 12.11.2001
% - added Langley's version of the quadratic transfer functions.
% revised pab 09.01.2001
% -simplified calculation of quadratic transfer functions when h=inf
% revised pab 09.01.2001
% - changed kurtosis so that kurtosis correspond to what kurt
% measures, i.e., kurtosis-3=0 for a Gaussian process
% by pab 01.03.2000


%error(nargchk(1,3,nargin))
narginchk(1,3)
% default options
opts.disp = 0;

if nargin<2||isempty(h),  h = S.h; end
if nargin<3||isempty(method), method = 'approximate'; end

S = spec2spec(S,'freq');
S = ttspec(S,'w');


S1 = S.S(:);
w = S.w(:);
m0 = trapz(w,S1);
Nw = length(w);

[Hs, Hd,Hdii] = qtf(S.w(:),h);

%return
%skew=6/sqrt(m0)^3*simpson(S.w,simpson(S.w,(Hs+Hd).*S1(:,ones(1,Nw))).*S1.');
Hspd = trapz(w',trapz(w,(Hs+Hd).*S1(:,ones(1,Nw))).*S1.');
switch lower(method(1))
  case 'a', %approx : Marthinsen, T. and Winterstein, S.R (1992) method
    if nargout>2
      ma = 2*trapz(S.w,Hdii.*S1);
    end
    sa   = sqrt(m0);
    skew = 6/sa^3*Hspd;
    kurt = (4*skew/3).^2+3;
  
  case 'q', % quasi method
    dw = diff(S.w(1:2));
    tmp1 =sqrt(S1(:,ones(1,Nw)).*(S1(:,ones(1,Nw)).'))*dw; 
    Hd = Hd.*tmp1;
    Hs = Hs.*tmp1;
    k = 6;
    stop = 0;
    while (~stop)
      E = eigs([Hd,Hs;Hs,Hd],[],k);
      %stop = (length(find(abs(E)<1e-4))>0 | k>1200);
      %stop = (any(abs(E(:))<1e-4) | k>1200);
      stop = (any(abs(E(:))<1e-4) | k>=min(2*Nw,1200));
      k = min(2*k,2*Nw);
    end
  
  
    m02=2*sum(E.^2); % variance of 2'nd order contribution 
  
    %Hstd = 16*trapz(S.w,(Hdii.*S1).^2);
    %Hstd = trapz(S.w,trapz(S.w,((Hs+Hd)+ 2*Hs.*Hd).*S1(:,ones(1,Nw))).*S1.');
    ma   = 2*trapz(S.w,Hdii.*S1);
    %m02  = Hstd-ma^2% variance of second order part
    sa   = sqrt(m0+m02);
    skew = 6/sa^3*Hspd;
    kurt = (4*skew/3).^2+3;
  case 'e', % Kac and Siegert eigenvalue analysis
    dw = diff(S.w(1:2));
    tmp1 =sqrt(S1(:,ones(1,Nw)).*(S1(:,ones(1,Nw)).'))*dw; 
    Hd = Hd.*tmp1;
    Hs = Hs.*tmp1;
    k = 6;
    stop = 0;
    %E2 = 1;
 
    while (~stop)
      [V,D] = eigs([Hd,Hs;Hs,Hd],[],k);
      E = diag(D);
      %stop = (length(find(abs(E)<1e-4))>0 | k>=min(2*Nw,1200));
      stop = (any(abs(E(:))<1e-4) | k>=min(2*Nw,1200));
      k = min(2*k,2*Nw);
    end
  
    
    h1 = sqrt(S1*dw/2);
    C  = (ctranspose(V)*[h1;h1]);
    
    E2 = E.^2;
    C2 = C.^2;
  
    ma   = sum(E);                     % mean 
    sa   = sqrt(sum(C2)+2*sum(E2));    % standard deviation
    skew = sum((6*C2+8*E2).*E)/sa^3;   % skewness
    kurt = 3+48*sum((C2+E2).*E2)/sa^4; % kurtosis
  otherwise, error('Method is not available')
end



return

function [Hs, Hd,Hdii]=qtf(w,h)
% QTF Quadratic Transfer Function
%
% CALL: [Hs, Hd, Hdii]=qtf(w,h)
%
%  Hs   = sum frequency effects
%  Hd   = difference frequency effects
%  Hdii = diagonal of Hd
%  w    = angular frequencies
%  h    = water depth


Nw = length(w);
g  = gravity;
kw = w2k(w,0,h,g);
[k1, k2] = meshgrid(kw);




if h==inf,% go here for faster calculations
  Hs   = 0.25*(abs(k1)+abs(k2));
  Hd   = -0.25*abs(abs(k1)-abs(k2));
  Hdii = zeros(size(w));
  return
% else
%   Hd = zeros(size(k1));
end

[w1, w2]=meshgrid(w);

msgId = 'MATLAB:divideByZero';
warning('off',msgId) %off % suppress warnings on division by zero

w12  = (w1.*w2);
w1p2 = (w1+w2);
w1m2 = (w1-w2);
k12  = (k1.*k2);
k1p2 = (k1+k2);
k1m2 = abs(k1-k2);
if 0, % Langley
  p1 = (-2*w1p2.*(k12*g^2-w12.^2)+...
      w1.*(w2.^4-g^2*k2.^2)+w2.*(w1.^4-g^2*k1.^2))./(4.*w12);
  p2= w1p2.^2.*cosh((k1p2).*h)-g*(k1p2).*sinh((k1p2).*h);
  
  Hs = -p1./p2.*w1p2.*cosh((k1p2).*h)/g-...
      (k12*g^2-w12.^2)./(4*g*w12)+(w1.^2+w2.^2)/(4*g);
  
  p3 = (-2*w1m2.*(k12*g^2+w12.^2)-...
      w1.*(w2.^4-g^2*k2.^2)+w2.*(w1.^4-g^2*k1.^2))./(4.*w12);
  p4= w1m2.^2.*cosh(k1m2.*h)-g*(k1m2).*sinh((k1m2).*h);
  
   
  Hd = -p3./p4.*(w1m2).*cosh((k1m2).*h)/g-...
      (k12*g^2+w12.^2)./(4*g*w12)+(w1.^2+w2.^2)/(4*g);  

else  % Marthinsen & Winterstein
  tmp1 = 0.5*g*k12./w12;
  tmp2 = 0.25/g*(w1.^2+w2.^2+w12);
  Hs   = (tmp1-tmp2+0.25*g*(w1.*k2.^2+w2.*k1.^2)./(w12.*(w1p2))).....
      ./(1-g*(k1p2)./(w1p2).^2.*tanh((k1p2).*h))+tmp2-0.5*tmp1; % OK
  
  tmp2 = 0.25/g*(w1.^2+w2.^2-w12); %OK
  Hd   = (tmp1-tmp2-0.25*g*(w1.*k2.^2-w2.*k1.^2)./(w12.*(w1m2))).....
      ./(1-g*(k1m2)./(w1m2).^2.*tanh((k1m2).*h))+tmp2-0.5*tmp1; % OK
end  

%tmp1 = 0.5*g*kw./(w.*sqrt(g*h));
%tmp2 = 0.25*w.^2/g;


Cg   = 0.5*g*(tanh(kw*h) +kw*h.*(1- tanh(kw*h).^2))./w; %Wave group velocity
Hdii = 0.5*(0.5*g*(kw./w).^2-0.5*w.^2/g+g*kw./(w.*Cg))....
      ./(1-g*h./Cg.^2)-0.5*kw./sinh(2*kw*h); % OK
Hd(1:Nw+1:end) = Hdii;

%k    = find(w1==w2);
%Hd(k) = Hdii;

% The NaN's occur due to division by zero. => Set the isnans to zero
Hdii(isnan(Hdii))=0;
Hd(isnan(Hd))=0;
Hs(isnan(Hs))=0;

warning('on',msgId)
return


