function [a,b]=fourier(t,x,T,M,N)
%FOURIER Returns Fourier coefficients.
%
%  CALL:  [a,b] = fourier(t,x,T,M);
%
%    a,b  = Fourier coefficients size M x P
%     t   = vector with N values indexed from 1 to N.
%     x   = vector or matrix of column vectors with data points size N x P.
%     T   = primitive period of signal, i.e., smallest period. 
%           (default T = diff(t([1,end]))
%     M-1 = no of harmonics desired (default M = N)
%     N   = no of data points (default length(t))
%
% FOURIER finds the coefficients for a Fourier series representation
% of the signal x(t) (given in digital form).  It is assumed the signal
% is periodic over T.  N is the number of data points, and M-1 is the
% number of coefficients.
%
% The signal can be estimated by using M-1 harmonics by:
%                    M
% x(i) = 0.5*a(1) + sum [a(n)*c(n,i) + b(n)*s(n,i)]
%                   n=2
% where
%   c(n,i) = cos(2*pi*(n-1)*t(i)/T)
%   s(n,i) = sin(2*pi*(n-1)*t(i)/T)
%
% Note that a(1) is the "dc value".
% Remaining values are a(2), a(3), ... , a(M).
% 
% Example
%  T = 2*pi;M=5;
%  t = linspace(0,4*T).'; x = sin(t);
%  [a,b] = fourier(t,x,T,M)
%
% See also  fft


%History:
% Revised pab 22.06.2001
%  -updated help header to wafo style.
%  -vectorized for loop.
%  -added default values for N,M, and T
%
% ME 244L Dynamic Systems and Controls Laboratory
% R.G. Longoria 9/1998
%

t = t(:);
if nargin<5||isempty(N), N = length(t);end
if nargin<4||isempty(M), M = N;        end
if nargin<3||isempty(T), T = diff(t([1,end]));end
szx = size(x);
P=1;
if prod(szx)==N,
  x = x(:);
elseif szx(1)==N
  P = prod(szx(2:end));
else
  error('Wrong size of x')
end


switch 0
  case -1,
       % Define the vectors for computing the Fourier coefficients
  %
  a = zeros(M,P);
  b = zeros(M,P);
  a(1,:) = simpson(x);

  %
  % Compute M-1 more coefficients
  tmp  = 2*pi*t(:,ones(1,P))/T;
  % tmp =  2*pi*(0:N-1).'/(N-1); 
  for n1 = 1:M-1,
    n = n1+1;
    a(n,:) = simpson(x.*cos(n1*tmp));
    b(n,:) = simpson(x.*sin(n1*tmp));
  end
  
  a = 2*a/N;
  b = 2*b/N;
  
  case 0,
     %
  a = zeros(M,P);
  b = zeros(M,P);
  a(1,:) = trapz(t,x);

  %
  % Compute M-1 more coefficients
  tmp  = 2*pi*t(:,ones(1,P))/T;
  % tmp =  2*pi*(0:N-1).'/(N-1); 
  for n1 = 1:M-1,
    n = n1+1;
    a(n,:) = trapz(t,x.*cos(n1*tmp));
    b(n,:) = trapz(t,x.*sin(n1*tmp));
  end
  a = a/pi;
  b = b/pi;
  
  case 1,
  % Define the vectors for computing the Fourier coefficients
  %
  a = zeros(M,P);
  b = zeros(M,P);
  
  %
  % Compute the dc-level (the a(0) component).
  %
  % Note: the index has to begin with "1".
  %
  
  a(1,:) = sum(x);

  %
  % Compute M-1 more coefficients
  tmp  = 2*pi*t(:,ones(1,P))/T;
  % tmp =  2*pi*(0:N-1).'/(N-1); 
  for n1 = 1:M-1,
    n = n1+1;
    a(n,:) = sum(x.*cos(n1*tmp));
    b(n,:) = sum(x.*sin(n1*tmp));
  end
  a = 2*a/N;
  b = 2*b/N;
case 2,
   % Define the vectors for computing the Fourier coefficients
  %
  a = zeros(M,P);
  b = zeros(M,P);
  a(1,:) = trapz(x);

  %
  % Compute M-1 more coefficients
  tmp  = 2*pi*t(:,ones(1,P))/T;
  % tmp =  2*pi*(0:N-1).'/(N-1); 
  for n1 = 1:M-1,
    n = n1+1;
    a(n,:) = trapz(x.*cos(n1*tmp));
    b(n,:) = trapz(x.*sin(n1*tmp));
  end
  
  a = 2*a/N;
  b = 2*b/N;
case 3
  % Alternative:  faster for large M, but gives different results than above.
  nper = diff(t([1 end]))/T; %No of periods given
  if nper == round(nper), 
    N1 = N/nper;
  else
    N1 = N;
  end

  % Fourier coefficients by fft
  Fcof1 = 2*ifft(x(1:N1,:),[],1);
  Pcor = [1; exp(sqrt(-1)*(1:M-1).'*t(1))]; % correction term to get
                                              % the correct integration limits
  Fcof = Fcof1(1:M,:).*Pcor(:,ones(1,P));
  a = real(Fcof(1:M,:));
  b = imag(Fcof(1:M,:));
end
return

%Old call: kept just in case

%
% Compute the dc-level (the a(0) component).
%
% Note: the index has to begin with "1".
%8
% a(1) = 2*sum(x)/N;
% 
% %
% % Compute M-1 more coefficients
% for n = 2:M,
%   sumcos=0.0;
%   sumsin=0.0;
%   for i=1:N,
%     sumcos = sumcos + x(i)*cos(2*(n-1)*pi*t(i)/T);
%     sumsin = sumsin + x(i)*sin(2*(n-1)*pi*t(i)/T);
%   end
%   a(n) = 2*sumcos/N;
%   b(n) = 2*sumsin/N;
% end
% 
% 
% return
% 
% 
