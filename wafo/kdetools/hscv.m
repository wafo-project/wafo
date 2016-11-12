function [h,hvec,score]=hscv(A,kernel,hvec)
%HSCV Smoothed cross-validation estimate of smoothing parameter.
%
% CALL: [hs,hvec,score] = hscv(data,kernel,hvec); 
% 
%   hs     = smoothing parameter
%   hvec   = vector defining possible values of hs
%             (default linspace(0.25*h0,h0,100), h0=hos(data,kernel))
%   score  = score vector
%   data   = data vector
%   kernel = 'gaussian'      - Gaussian kernel the only supported
%                               
%  Note that only the first 4 letters of the kernel name is needed.
%  
% Example: 
%   data = rndnorm(0,1,20,1)
%   [hs hvec score] = hscv(data,'epan');
%   plot(hvec,score) 
%
%   close all;
%
% See also  hste, hbcv, hboot, hos, hldpi, hlscv, hstt, kde, kdefun


% tested on : matlab 5.2
% history:
% revised pab Aug 2005
% -added support for other kernels by scaling 
% (this is an ad hoc solution, which may be improved in future)
% revised pab dec2003  
% revised pab 20.10.1999
%   updated to matlab 5.2
% changed input arguments
% improoved gridding
% taken from kdetools     Christian C. Beardah 1995 

%  Wand,M.P. and Jones, M.C. (1986) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 75--79

% TODO % Add support for other kernels than Gaussian  
%A=A(:);
[n, d] = size(A);
if (n==1) && (d>1),
  A=A.';
  n=d;
  d=1;
end

if nargin<2||isempty(kernel),
  kernel='gauss';
end;

% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
[mu2,R] = kernelstats(kernel);
AMISEconstant = (8 * sqrt(pi) * R / (3 * mu2 ^ 2 * n)) ^ (1. / 5);
STEconstant = R /(mu2^(2)*n);

sigmaA = hns(A, kernel)/AMISEconstant;
if nargin<3||isempty(hvec),
  %H=hos(A,kernel);
  H = AMISEconstant / 0.93;
  hvec=linspace(0.25*H,H,100);
else
  hvec=abs(hvec);
end;
  
steps=length(hvec);
score = zeros(1,steps);

inc  = 128;
nfft = inc*2;

xmin=min(A);    % Find the minimum value of A.
xmax=max(A);    % Find the maximum value of A.
xrange=xmax-xmin; % Find the range of A.

% sigmaA = std(A);
% iqr = abs(diff(qlevels2(A,[75 25]),1,1+(d==1))); % interquartile range
% k = find(iqr>0);
% if any(k)
%   sigmaA(k) = min(sigmaA(k), iqr(k)/1.349);
% end



% xa holds the x 'axis' vector, defining a grid of x values where 
% the k.d. function will be evaluated and plotted.

ax1=xmin-xrange/8;
bx1=xmax+xrange/8;

kernel2 = 'gauss';
[mu2,R] = kernelstats(kernel2);
STEconstant2 =  R /(mu2^(2)*n);

h = zeros(1,d);
hvec = hvec*(STEconstant2/STEconstant)^(1/5);
for dim = 1:d
  s  = sigmaA(dim);
  ax = ax1(dim)/s;
  bx = bx1(dim)/s;
  xa = linspace(ax,bx,inc).'; 
  xn = linspace(0,bx-ax,inc);
  datan = A(:,dim)/s;
  c = gridcount(datan,xa);

  %deltax = (bx-ax)/(inc-1);
  %xn     = (0:(inc-1))*deltax;
  %binx=floor((A(:,dim)-ax)/deltax)+1;
  % Obtain grid counts
  %c = full(sparse(binx,1,(xa(binx+1)-A(:,dim)),inc,1)+...
  %	   sparse(binx+1,1,(A(:,dim)-xa(binx)),inc,1))/deltax;

  [k40,k60,k80,k100] = deriv(0,kernel2);
  psi8  = 105/(32*sqrt(pi)); %*s^9);
  psi12 = 3465/(512*sqrt(pi)); %*s^13);
  g1 = (-2*k60/(mu2*psi8*n))^(1/9);
  %g1 = sqrt(2)*s*(2/(7*n))^(1/9)
  g2 = (-2*k100/(mu2*psi12*n))^(1/13);
  %g2 = sqrt(2)*s*(2/(11*n))^(1/13)

  [kw4,kw6] = deriv(xn/g1,kernel2);% Obtain the kernel weights.
  kw = [kw6,0,kw6(inc:-1:2)].';% Apply 'fftshift' to kw.
  z  = real(ifft(fft(c,nfft).*fft(kw)));% Perform the convolution.

  psi6 = sum(c.*z(1:inc))/(n^2*g1^7);

  % Obtain the kernel weights.

  [kw4,kw6,kw8,kw10]=deriv(xn/g2,kernel2);
  kw=[kw10,0,kw10(inc:-1:2)]';% Apply 'fftshift' to kw.
  z=real(ifft(fft(c,nfft).*fft(kw)));% Perform the convolution.

  psi10=sum(c.*z(1:inc))/(n^2*g2^11);
  g3 = (-2*k40/(mu2*psi6*n))^(1/7);
 % g3 = (-6/(sqrt(2*pi)*psi6*n))^(1/7)
  g4 = (-2*k80/(mu2*psi10*n))^(1/11);
  %g4=(-210/(sqrt(2*pi)*psi10*n))^(1/11)
  
  % Obtain the kernel weights.
  
  kw4=deriv(xn/g3,kernel2);
  
  % Apply 'fftshift' to kw.
  
  kw=[kw4,0,kw4((inc:-1:2))]';
  
  % Perform the convolution.
  
  z=real(ifft(fft(c,nfft).*fft(kw)));
  
  psi4=sum(c.*z(1:inc))/(n^2*g3^5);
  
  % Obtain the kernel weights.
  
  [kw4,kw6,kw8]=deriv(xn/g4,kernel2);
  
  % Apply 'fftshift' to kw.
  
  kw=[kw8,0,kw8(inc:-1:2)]';
  
  % Perform the convolution.
  
  z=real(ifft(fft(c,nfft).*fft(kw)));
  
  psi8=sum(c.*z(1:inc))/(n^2*g4^9);

  C=(441/(64*pi))^(1/18)*(4*pi)^(-1/5)*psi4^(-2/5)*psi8^(-1/9);
  
  M=datan*ones(size(datan'));
  
  Y=(M-M');
  
  for i=1:steps,
    
    g=C*n^(-23/45)*hvec(i)^(-2);
    
    sig1=sqrt(2*hvec(i)^2+2*g^2);
    
    sig2=sqrt(hvec(i)^2+2*g^2);
    
    sig3=sqrt(2*g^2);
    
    term2=sum(sum(mkernel(Y/sig1,kernel2)/sig1-2*mkernel(Y/sig2,kernel2)/sig2+mkernel(Y/sig3,kernel2)/sig3));
    
    score(i)=1/(n*hvec(i)*2*sqrt(pi))+term2/n^2;
    
  end;

  [L,I]=min(score);
  
  h(dim) = sigmaA(dim)*hvec(I);
  
   % Kernel other than Gaussian scale bandwidth
  h(dim)  = h(dim)*(STEconstant/STEconstant2)^(1/5);
end

hvec = hvec*(STEconstant/STEconstant2)^(1/5);

return
data =[ 0.75355792,  0.72779194,  0.94149169,  0.07841119,  2.32291887,...
    1.10419995,  0.77055114,  0.60288273,  1.36883635,  1.74754326,...
    1.09547561,  1.01671133,  0.73211143,  0.61891719,  0.75903487,...        
    1.8919469 ,  0.72433808,  1.92973094,  0.44749838,  1.36508452]