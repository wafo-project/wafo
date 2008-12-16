function [x, xder]=cov2sdat(R,np,iseed)
%COV2SDAT Simulates a Gaussian process and its derivative
%         from covariance func.
%           
%   CALL: [xs, xsder ] = cov2sdat(R,[np cases],iseed);
% 
%         xs    = a cases+1 column matrix  ( t,X1(t) X2(t) ...).
%         xsder = a cases+1 column matrix  ( t,X1'(t) X2'(t) ...). 
%         R     = a covariance structure
%         np    = giving np load points.  (default length(S/R)-1=n-1).
%                 If np>n-1 it is assummed that R(k)=0 for all k>n-1
%         cases = number of cases, i.e. number of replicates (default=1) 
%         iseed = starting seed number for the random number generator 
%                 (default none is set)
%
%
%  COV2SDAT performs a Fast and exact simulation of stationary zero mean 
%  Gaussian process through circulant embedding of the Covariance matrix.
%  
%  If the covariance structure  has a non-empty field .tr, then the
%  transformation is  applied to the simulated data, the result is a 
%  simulation of a transformed Gaussian process.
%
% Example:
%  np = 100; 
% [x1, x2] = cov2sdat(spec2cov(jonswap),np);
% waveplot(x1,'r',x2,'g',1,1)  
%
% See also  spec2sdat, tranproc

% Reference 
% C.R Dietrich and G. N. Newsam (1997)
% "Fast and exact simulation of  stationary 
%  Gaussian process through circulant embedding 
%  of the Covariance matrix" 
%  SIAM J. SCI. COMPUT. Vol 18, No 4, pp. 1088-1107


% tested on: Matlab 5.3
% History:
% Revised pab Feb2004
% - changed seed to state   
% revised pab 12.10.2001
% - added example
% - the derivative, xder, is now calculated correctly using fft.
%   derivate.m is no longer needed!
% revised pab 11.10.2001
% - added call to lagtype.m
% - added a new solution on removing high frequency noise due to the fact
%  that the fix of 16.01.2000 was not sufficient.
% revised by es 24.05.00 Some changes in help  
% revised pab 13.03.2000
% -commented out warning messages
% revised pab 24.01.2000
% -added iseed
% revised pab 16.01.2000
%  -added a fix in line 134: truncating the values of the spectral density to
%     zero if S/max(S) > trunc. This is to ensure that high frequency 
%     noise is not added to the simulated timeseries.
%  -fixed transformation for the derivatives
% revised pab 12.11.1999
%  -fixed transformation
% revised pab 12.10.1999
%    last modified by Per A. Brodtkorb 19.08.98

% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates
nugget=0;%10^-12;

ACF=R.R;
[n,m]=size(ACF);

if n<m
 b=m;m=n;n=b; 
 ACF=ACF';
end

if n<2, 
  error('The input vector must have more than 2 elements!')
end

%istime=1;
switch m
 case 1, % OK
 otherwise, error('Only capable of 1D simulation')          
end

[y I]=max(ACF);
switch I,
  case 1,% ACF starts with zero lag
  otherwise, error('ACF does not have a maximum at zero lag')  
end
% new call pab 11.10.2001
ltype = lagtype(R);
%names=fieldnames(R);
%ind=find(strcmp(names,'x')+strcmp(names,'y')+strcmp(names,'t'));
      %options are 'f' and 'w' and 'k'      
%ltype=lower(names{ind});

ti=R.(ltype);
if isempty(ti)
  dT=1;
else
  dT=ti(2)-ti(1); % ti 
end

if nargin<2||isempty(np)
  np=n-1;
  cases=1;
else
  switch  length(np) 
    case 1, cases=1; 
    case 2, cases=np(2); np=np(1);
    otherwise, error('Wrong input. Too many arguments')
  end
end
if nargin<3||isempty(iseed)
else
  try
    randn('state',iseed);
  catch
    randn('seed',iseed); % set the the seed  
  end
end
x=zeros(np,cases+1);


if nargout==2
  xder=x;
end
  
% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates
ACF(1)=ACF(1)+nugget;
   
% Fast and exact simulation of simulation of stationary
% Gaussian process throug circulant embedding of the
% Covariance matrix
if (abs(ACF(n))>eps),% assuming ACF(n+1)==0
  m2=2*n-1;
  nfft=2^nextpow2(max(m2,2*np));
  ACF=[ACF(1:n);zeros(nfft-m2,1);ACF(n:-1:2)];
  if 0, %(n<np),
    disp('Warning: I am now assuming that ACF(k)=0 ')
    disp('for k>MAXLAG.')       
  end
else % ACF(n)==0
  m2=2*n-2;
  nfft=2^nextpow2(max(m2,2*np));
  ACF=[ACF(1:n);zeros(nfft-m2,1);ACF(n-1:-1:2)];
end
%m2=2*n-2;

S=real(fft(ACF,nfft));% periodogram


[maxS I]=max(S(:));
k=find(S<0);
if any(k),%
  %  disp('Warning: Not able to construct a nonnegative circulant ')
  %  disp('vector from the ACF. Apply the parzen windowfunction ') 
  %  disp('to the ACF in order to avoid this.')
  %  disp('The returned result is now only an approximation.')
  
  % truncating negative values to zero to ensure that 
  % that this noise is not added to the simulated timeseries
 
  S(k)=0; 
  % New call pab 11.10.2001	    
  ix = min(k((k>2*I)));
  if any(ix),
    % truncating all oscillating values above 2 times the peak 
    % frequency to zero to ensure that 
    % that high frequency noise is not added to 
    % the simulated timeseries.
    S(ix:end-ix) =0;
  end
end


trunc=1e-5;
k=find(S(I:end-I)/maxS<trunc);
if any(k),%
 S(k+I-1)=0; 
 % truncating small values to zero to ensure that 
 % that high frequency noise is not added to 
 % the simulated timeseries	        
end

cases2=ceil(cases/2);
% Generate standard normal random numbers for the simulations
% -----------------------------------------------------------
if 1,
  epsi=randn(nfft,cases2)+sqrt(-1)*randn(nfft,cases2);
  Ssqr=sqrt(S/(nfft)); %sqrt(S(wn)*dw )
  ephat=epsi.*Ssqr(:,ones(1,cases2));
  y=fft(ephat,nfft);
  x(:,2:cases+1)=[real(y(3:np+2,1:cases2)) imag(y(3:np+2,1:floor(cases/2)))]; 
else
  epsi=randn(nfft/2+1,cases2)+sqrt(-1)*randn(nfft/2+1,cases2);
  epsi(1,:)=real(epsi(1,:));
  epsi(end,:)=real(epsi(end,:));
  epsi=[epsi; conj(epsi((end-1):-1:2,:))];  
end


x(:,1)=(0:dT:(dT*(np-1)))';

if nargout==2,
  
    % xder(:,1:(cases+1))=derivate((-2*dT:dT:(dT*(np+1)))', ...
    %	[real(y(1:np+4,1:cases2)) imag(y(1:np+4,1:floor(cases/2)))]); 
  
    % New call: pab 12.10.2001
    Ssqr  = Ssqr.*[(0:(nfft/2))';-((nfft/2-1):-1:1)']*2*pi/nfft/dT;  
    ephat = epsi.*Ssqr(:,ones(1,cases2));
    y     = fft(ephat,nfft);    
    xder(:,2:(cases+1))=[imag(y(3:np+2,1:cases2)) -real(y(3:np+2,1:floor(cases/2)))];      
    xder(:,1)=x(:,1);
  
end

if isfield(R,'tr') && ~isempty(R.tr)
  disp('   Transforming data.')
  g=R.tr;
  G=fliplr(g); % the invers of g
  if nargout==2, % gaus2dat
    for ix=1:cases
      tmp=tranproc([x(:,ix+1) xder(:,ix+1)],G); 
      x(:,ix+1)=tmp(:,1);
      xder(:,ix+1)=tmp(:,2);
    end
  else
    for ix=1:cases % gaus2dat
      x(:,ix+1)=tranproc(x(:,ix+1),G); 
    end
  end
end
