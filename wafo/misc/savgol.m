function c = savgol(nl,nr,m,ld,np)
%SAVGOL  Savitzky-Golay filter coefficients.
%
% CALL:  c = savgol(nl,nr,m,ld,np);
%
%  c     = vector with Savitzky-Golay filter coefficients.
%  nl,nr = number of leftward (past) and rightward (future) data points
%          used, respectively, making the total number of data points
%          used nl+nr+1. (nl+nr>=m)
%  m     = order of the smoothing polynomial, also equal to the highest
%          conserved moment; usual values are m = 2 or m = 4.
%  ld    = order of the derivative desired. (default ld=0 for smoothed
%          function) (ld<=m)
%  np    = length of c  (np>=nr+nl+1) (defalt 2*max(nr,nl)+1) 
% 
% The idea of Savitzky-Golay filtering is to find filter coefficients Cn
% that preserves higher moments, i.e., to approximate the underlying
% function within the moving window not by a constant (whose estimate is
% the average), but by a polynomial of higher order, typically quadratic
% or quartic. For each point Y(i), we least-squares fit a polynomial of
% order M to all NL+NR+1 points in the moving window and then set YF(i)
% to the value of this polynomial at position I, i.e., 
%
%             nr
%    YF(i) = sum Cn(n)*Y(i+n)
%           n=-nl
% 
% SAVGOL Returns the filter coefficients C(1:np), in wrap-around order,
% i.e., C(1) = Cn(0), C(2)=Cn(-1),...C(Nl+1)=Cn(-Nl), and C(Np) = Cn(1),
% C(Np-1) = Cn(2),...C(Np-Nr) = Cn(Nr), which is consistent for use with
% fft to perform the convolution.
%
% Note the restrictions:  np >= nr+nl+1 > m >= ld
%
% Example:
%   x = linspace(0,1);
%   y = exp(x)+1e-1*randn(size(x));
%   nl=2; nr=2; m=2; ld=0;
%   c = savgol(nl,nr,m,ld,length(x));
%   yf = real(ifft(fft(y).*fft(c))); % convolution with pronounced end effects
%   yf2 = convlv(y,c);               % convolution with less end effects
%   plot(x,y,x,yf,x,yf2);
%   ld =1; m =4;
%   c = savgol(nl,nr,m,ld);          % differentiation filter
%   dyf = convlv(y,c)*gamma(ld+1)/(x(2)-x(1))^ld; % Derivative
%   ix = nl+1:length(x)-nr;                       % for all these ix
%   semilogy(x(ix),abs(y(ix)-dyf(ix))); % Error of the derivative
%
%   close all;
%
% See also  cssmooth

% Reference 
% William H. Press, Saul Teukolsky, 
% William T. Wetterling and Brian P. Flannery (1997)
% "Numerical recipes in Fortran 77", Vol. 1, pp 644-649

% History
% by pab 2000


error(nargchk(3,5,nargin))
if nargin<5||isempty(np),np=2*max(nl,nr)+1;end
if nargin<4||isempty(ld),ld=0;end

if (np<nl+nr+1)
  error('WAFO:SAVGOL','np must be larger than nl+nr = %d or nl+nr less than np = %d',nr+nl,np)
end
if nl<0, 
  warning('WAFO:SAVGOL','nl must be positive or zero')
  nl = max(nl,0);
end
if nr<0, 
  warning('WAFO:SAVGOL','nr must be positive or zero')
  nr = max(nr,0);
end
if ld>m 
  error('WAFO:SAVGOL','m must be larger than ld = %d ',ld)
end
if nl+nr<m,
  error('WAFO:SAVGOL','m must be smaller than nl+nr = %d',nl+nr)
end

asiz = [m+1,m+1];
a    = zeros(asiz);
for ipj=0:2*m  
  %Set up the normal equations of the desired least-squares fit.   
  tmp = sum((1:nr).^ipj) + (ipj==0) + sum((-1:-1:-nl).^ipj); 
  
  mm = min(ipj,2*m-ipj);
  imj=-mm:2:mm;
  ind = sub2ind(asiz, 1+(ipj+imj)/2,1+(ipj-imj)/2);
  a(ind)=tmp;
  %for imj=-mm:2:mm 
  %  a(1+(ipj+imj)/2,1+(ipj-imj)/2)=tmp;
  %end 
end 

% Solve them by LU decomposition.
[L,U] = lu(a);

% Right-hand side vector is unit vector, 
% depending on which derivative we want.
b       = zeros(m+1,1);
b(ld+1) = 1;

% Backsubstitute, giving one row of the inverse matrix.

d = (U\(L\b)).';


%Zero the output array (it may be bigger than number of coefficients). 
c = zeros(1,np);

if 0,
  for k=-nl:nr 
    %Each Savitzky-Golay coefficient is the dot product
    %of powers of an integer with the inverse matrix row.
    tmp = sum(d(1:m+1).*[1 k.^(1:m)]);
    kk = mod(np-k,np)+1; %Store in wrap-around order.
    c(kk) = tmp;
  end 
else
  %Each Savitzky-Golay coefficient is the dot product
  %of powers of an integer with the inverse matrix row.
  k    = (-nl:nr ).';
  Nk   = length(k);
  tmp0 = repmat(k,1,m+1).^repmat(0:m,Nk,1);
  tmp  = sum(d(ones(Nk,1),1:m+1).*tmp0,2)';
  kk   = mod(np-k,np)+1; %Store in wrap-around order.
  c(kk) = tmp;
end
return

