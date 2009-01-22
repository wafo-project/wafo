function [corr,c0] = dat2cor(xn,varargin)
%DAT2COR Estimate auto correlation function from data.
%
% CALL:  r = dat2cor(x,L,plotflag,dT,flag)
%
%         r = Covariance structure containing:
%                 R = ACF vector length L+1  
%                 t = time lags  length L+1  
%             stdev = estimated large lag standard deviation of the estimate
%                     assuming x is a Gaussian process:
%                     if r(k)=0 for all lags k>q then an approximation 
%                     of the variance for large samples due to Bartlett
%                     var(corr(k))=1/N*(1+2*r(1)^2+2*r(2)^2+ ..+2*r(q)^2)
%                     for  k>q and where  N=length(x). Special case is
%                     white noise where it equals 1/N for k>0
%                 h = water depth (default inf)
%                tr = [], transformation 
%              type = 'none'
%              norm = 1 indicating that R is normalized
%
%        x  = a vector or a two column data matrix with 
%             sampled times and values. (length n)
%        L  = the maximum time-lag for which the ACF is estimated.
%	           (Default L=n-1)        
%  plotflag = 1 then the ACF is plotted vs lag
%	           2 then the ACF is plotted vs lag in seconds
%	           3 then the ACF is plotted vs lag and vs lag (sec
%        dT = time-step between data points (default xn(2,1)-xn(1,1) or 1 Hz).
%      flag = 'biased'  : scales the raw correlation by 1/n. (default)
%             'unbiased': scales the raw correlation by 1/(n-abs(k)), 
%                         where k is the lag.
%
%  Note: - flag may be given anywhere after x.
%        - x may contain NaN's (i.e. missing values).
%
% Example: 
%   x = load('sea.dat');
%   rf = dat2cor(x,150,2)
%
% See also  dat2cov, covplot

       
%(If you do have the signal toolbox you may modify this file)

% Tested on: Matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised jr 02.04.2001
%  - added example, updating help 
% revised pab 09.10.2000
%  - added secret output c0 = variance
%  - made computations faster for the case Ncens < n
%  - added flag option, dat2corchk
%  -fixed a bug: forgot to subtract the mean from x
% revised pab 12.08.99
% - changed code to handle covariance structure/class
% modified by Per A. Brodtkorb 28.10.98 
% enabled to handle missing data 


[x,L,plotflag,dT,flag,n] = dat2corchk(varargin,xn,nargout);

corr       = createcov(0,'t');
corr.R     = zeros(L+1,1);
corr.t     = zeros(L+1,1);
corr.stdev = zeros(L+1,1);
corr.h     = inf;
corr.norm  = 1; % normalized

%----------------------------------------------
%  
%  Signal toolbox users may modify below
%
%----------------------------------------------

indnan = isnan(x);
if any(indnan)
   x     = x-mean(x(~indnan)); % remove the mean pab 09.10.2000
   indnan = find(indnan);
   Ncens = n-length(indnan);
else
   indnan = [];
   Ncens  = n;
   x      = x-mean(x);
end

if 0 %Ncens<n, %Old call for missing data (slow)   
   r=zeros(L+1,1);
   % unbiased result, i.e. divide by n-abs(lag)
   for tau=0:L,
      Ntau=n-tau;
      ix=1:Ntau; % constructing a vector if indices
      %vector multiplication is faster than a for loop in MATLAB
      tmp=x(ix).*x(ix+tau);
      if strcmpi(flag,'unbiased'),
         r(tau+1)=mean(tmp(~isnan(tmp))); % unbiased estimate
      else
         r(tau+1)=sum(tmp(~isnan(tmp)))/Ncens; % biased estimate
      end
   end
else 
   if any(indnan)
      x(indnan)=0; % pab 09.10.2000 much faster for censored samples
   end
   if 1, % If you haven't got the signal toolbox use this
     nfft=2^nextpow2(n) ;
     Rper=abs(fft(x,nfft)).^2/Ncens; % Raw periodogram
     %plot(Rper)%,size(Rper),n,nfft
          
     r=real(fft(Rper))/nfft ; %ifft=fft/nfft since Rper is real!
      
     if strcmpi(flag,'unbiased'),
       % unbiased result, i.e. divide by n-abs(lag)
       r=r(1:L+1)*Ncens./ (Ncens:-1:Ncens-L)';
     %else  % biased result, i.e. divide by n
     %  r=r(1:L+1)*Ncens/Ncens;
     end

  else % If you have got the signal toolbox use this 
    r=xcorr(x,L,flag);
    r=(r((L+1):end))'; 
  end  
end

c0     = r(1);
corr.R = r(1:(L+1))/c0; %normalizing
corr.t = (0:L)'*dT;

corr.stdev=sqrt(1/Ncens*[ 0; 1 ;(1+2*cumsum(corr.R(2:L).^2))]);

if ((plotflag>0) ),
  covplot(corr,L,plotflag)
end

return

function [x,L,plotflag,dT,flag,n] = dat2corchk(P,xn,Na)
% DAT2CORCHK Helper function for dat2cor.
%
% CALL  [x2,L,plotflag,dT,flag,n]=dat2corchk(P,x1) 
%
%   P = the cell array P of input arguments (between 0 and 7 elements)
%   x = a vector.


x=xn;
[n m]= size(x);

if n<m
 b=m;m=n;n=b; 
 x=x';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

% Default values
%~~~~~~~~~~~~~~~~  
flag = 'biased';
plotflag = 0;
dT=[]; % sampling frequency unknown
L=n-1;
   % if n<400
%    L=floor(n/2);
%  else
%    L=floor(12*sqrt(n));
%  end




Np=length(P);
strix=zeros(1,Np);
for ix=1:Np, % finding symbol strings 
 strix(ix)=ischar(P{ix});
end


k=find(strix);
if any(k)
  flag = P{k(1)};
  Np   = Np-1;
  P={P{find(~strix)}};
end


if Np>0 && ~isempty(P{1}), L=P{1}; end
if Np>1 && ~isempty(P{2}), plotflag=P{2}; end
if Np>2 && ~isempty(P{3}), dT=P{3}; end

if (Na == 0)&&plotflag==0
  plotflag=3;
end
if plotflag>3, 
  error('Invalid option. Plotflag can only be 0,1,2 or 3')
end
 switch m
 case 1, x=x(:); if isempty(dT),dT=1;end
 case 2,  dT=x(2,1)-x(1,1);x=x(:,2);
 otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ') 
              
end
return
