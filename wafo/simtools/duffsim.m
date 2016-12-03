function X = duffsim(T,dt,z,a,b,alf,iseed,method)
% DUFFSIM Generates a sample path of a harmonic oscillator 
%
% CALL:  X = duffsim(T,dt,z,a,b,alf,iseed,method);
%
%        X   = a three column matrix with time in the first, the simulated 
%              process in the second and the derivative of the simulated 
%              process in the third column.
%        T   = the maximum time.
%        dt  = the time step.
%        z   = parameter of the equation for the oscillator.
%  a,b,alf   = parameters defining type of oscillator. 
%              (Default values 0,-1, and 2, respectively).
%     iseed  = seed to the random generator. (Default [] (random seed))
%     method = String defining random generator
%              'state'  : Period of 2^1492 (default)
%              'twister': Period of (2^19937-1)/2
%              'fortran': Fortran implementation (Fast, but period of random generator smaller)
%              
%    DUFFSIM generates a sample path of a harmonic oscillator with a nonlinear spring
%    driven by a white-noise process described by the equation:
% 
%        X''(t) + 2 z X'(t) + b X(t) + a X(t)^3 = 2 sqrt(z) W'(t),
%
%    where  W'(t) is the white-noise process and z,b,a  are constants. 
%    If ALF=2 then W'(t), is a Gaussian white noise.
%    If 0<ALF<2 then W'(t) is alpha-stable white noise . 
%
%    Important parameter values:
%    0<z<1, a=0, b=1 : Normalized linear oscillator, Var(X(t))=Var(X'(t))=1.
%    a=b=0, alf=2    : X'(t) is the Ornstein-Uhlenbeck process
%    a,z>0, b=-1     : Duffing oscillator.  
%    The simulation technique is Euler's discretization scheme.
%
% Example:
%  z = 0.7; a= 0; b = 1;alf = 2;
%  T = 10000; dt = 0.5; L = 120;
%  x = duffsim(T,dt,z,a,b,alf);
%  plot(x(:,1),x(:,2))
%  S   = dat2spec(x(:,1:2),L);
%  Hm0 = spec2char(S,'hm0');
%  w   = [0.01 2*pi/dt 275];
%  Se  = oscspec(w,z/sqrt(b),sqrt(b),2*sqrt(z/b));
%  plotspec(Se),hold on,
%  plotspec(S,'r'),hold off, shg
%
%  close all;
% 
% See also oscspec

% History: 
% -revised pab Feb 2007
%  -Added a matlab  version of the oscillator.
%  -Added iseed + method.
%  -Added example
% Adapted from WAT. 
% revised jr: 00.05.16
% - updated final loading
% - updated information

if ( (z<=0) || (z>1) )
  error('WAFO:DUFFSIM','   Parameter z not in (0,1).')
end
if nargin<8 || isempty(method)
 method = 'state'; % Period of 2^1492
 %method = 'twister'; % Period of (2^19937-1)/2
 % method = 'fortran'; % unknown period but smaller than above.
end
if nargin<7 || isempty(iseed)
 iseed=[];
elseif ~strncmpi(method,'fortran',1)
  try
    rand(method,iseed)
  catch
    warning('WAFO:DUFFSIM','Uknown random generator! Using multiplicative congruential algorithm with period of 2^31-2.')
    rand('seed',iseed); % set the the seed  
  end
end

if nargin<6 || isempty(alf)
   alf=2;
end
if nargin<5 || isempty(b)
   b=-1;
end

if nargin<4 || isempty(a)
   %b=1;
   a=0;
end
N=(floor(T/dt));
if (N>500000)
  error('WAFO:DUFFSIM','Number of points requested very large or time step  dt  is too small, break.')
end

if strncmpi(method,'fortran',1)
  try
    X = simduff_f77(N,dt,z,a,b,alf,iseed);
  catch
    warning('WAFO:DUFFSIM','Failed calling fortran version. Calling matlab version instead!')
    X = simduff(N,dt,z,a,b,alf);
  end
else
  X = simduff(N,dt,z,a,b,alf);
end

end % duffsim

function X = simduff_f77(N,dt,z,a,b,alf,iseed)
% SIMDUFF_F77 Harmonic oscillator Fortran version


if exist('simduff.in','file')
  delete simduff.in,
end
disp('   Writing data.')
data = [100*N 0.01*dt z a b alf];
if isempty(iseed)
  iseed = floor(1e8+rand*899999999);
end
fid = fopen('simduff.in','w');
fprintf(fid,'%6.0f %7.5f %7.5f %7.5f %7.5f %7.5f\n',data);
fprintf(fid,'%10.0f\n',iseed);
fclose(fid);
disp('   Starting Fortran executable.')
dos([wafoexepath 'simduff.exe']);

disp('   Loading data.')
X = load('out.dat');
delete simduff.in
end % simduff_f77

function X = simduff(N,dt,z,a,b,alf)
% SIMDUFF  Harmonic oscillator Matlab version

Nsim   =  N*100;
dt1    = dt/100;
X      = zeros(N,3);
X(:,1) = (1:N).'*dt;

s  = 2*sqrt(z)*exp(1./alf*log(dt1));
if alf>1.99999
  random = gaussrnd([Nsim+1,1]);
else
  random = stable([Nsim+1,1]);
end

%random = randgen(Bsz);
X1   = random(1);
if (a>0.01)
  X1 = X1/sqrt(a);
end
Xdot1 = random(2);
kp = 1;
h = fwaitbar(0,[],'this may take a while');
for i=2:Nsim
  X2    = X1+Xdot1*dt1;
  ix= i+1;
  Xdot2 = Xdot1+(-b*X1-a*X1.^3-2*z*Xdot1)*dt1 + s*random(ix);
  X1    = X2;
  Xdot1 = Xdot2;
  if (i/100==kp)
    X(kp,2:3) = [X2,Xdot2];
    kp = kp+1;
    if mod(kp,100)==0
      fwaitbar(kp/N,h)
    end
  end
end
close(h)
   
  
  function r = gaussrnd(sz)
    %GAUSSRND Random numbers from a Gaussian distribution
    u1 = rand(sz);
    u2 = rand(sz);
    r = sqrt(-2*log(u1)).*cos(2*pi*u2);
    
  end % gaussrnd

  function y = stable(sz)
    %data pi/3.141592654/
  if nargin<1
    sz = 1;
  end
  v = pi*(rand(sz)-0.5);
  xk = alf*v;
  w  = -log(rand(sz))./cos(v-xk);
  hh = w.*cos(v);
  y  = sin(xk).*exp(-log(hh)/alf).*w/sqrt(2);
  return
  end % stable
end % simduff

