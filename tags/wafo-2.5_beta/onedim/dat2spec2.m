function S = dat2spec2(xn,L,g,plotflag,wdef,dT,chopOffHighFreq)
%DAT2SPEC2 Estimate one-sided spectral density, version 2. 
%	  using a Parzen window function on 
%	  the estimated autocovariance function.
% 
% CALL:  S = dat2spec2(x,L,g,plotflag,wdef,dT)
%
%         S = A structure containing:
%             S    = spectral density
%             w    = angular frequency
%             tr   = transformation g
%             h    = water depth (default inf)
%             type = 'freq'
%             note = Memorandum string
%             date = Date and time of creation 
%             L    = maximum lag size of the window function. 
%             Bw   = Bandwidth of the smoothing window which is used 
%                    in the estimated spectrum. (rad/sec or Hz)
%            
%
%         x =  m column data matrix with sampled times in the first column
%              and values the next columns.    
%
%         L = maximum lag size of the window function. 
%             If no value is given, the lag size is set to
%             be the lag where the auto correlation is less than 
%             2 standard deviations. (maximum 300) 
%	                     
%         g = the transformation assuming that x is a sample of a 
%             transformed Gaussian process. If g is empty then
%             x  is a sample of a Gaussian process (Default)
%
%  plotflag = 1 plots the spectrum, S, 
%             2 plot 10log10(S) and
%             3 plots both the above plots
%
%        wdef = 'parzen'   then a Parzen window is used (default).
%               'hanning'  then a Hanning window is used. 
%
%         dT  = sampling interval. Default is dT=x(2,1)-x(1,1) or 1Hz.
%
%  As L decreases the estimate becomes smoother and Bw increases. 
%  If we want to resolve peaks in S which is Bf (Hz or rad/sec) apart 
%  then Bw < Bf. 
%  
% See also   dat2tr, dat2cov, dat2spec

% May chop off high frequencies in order to 
% get the same irregularity factor in the spectrum as 
% in x.


% References:
% Georg Lindgren and Holger Rootzen (1986)
% "Stationära stokastiska processer",  pp 173--176.
% 
% Gareth Janacek and Louise Swift (1993)
% "TIME SERIES forecasting, simulation, applications",
% pp 75--76 and 261--268
%
% Emanuel Parzen (1962),
% "Stochastic Processes", HOLDEN-DAY,
% pp 66--103


% Tested on: Matlab 5.3
% history:
% Modified by Per A. Brodtkorb 14.08.98,25.05.98
%  - add a nugget effect to ensure that round off errors
%    do not result in negative spectral estimates
% Modified by svi 29.09.99
%  - The program is not estimating transformation g any more.
% modified pab 22.10.1999
%  - fixed so that x in fact can be a m column matrix
%  - updated info put into the spectral structure.
%  - updated help header
% modified pab 03.11.1999
%  - fixed a bug: line 152 wrong array dim. when m>2
% modified jr 29.01.2000
%  - modified from dat2spec, new name: dat2spec2
%  - lines related to confidence bands removed
%  - call to psd removed
% Revised pab Dec 2003
% -fixed a bug for wdef and updated helpheader accordingly.
%  

nugget=10^-12;

freqtype='w'; %options are 'f' and 'w'

xx=xn;
istime=1;
[n m]= size(xx);

if min(m,n)==1, 
  xx=[ (1:n)' xx(:)];istime=0;
  n=max(m,n);
  m=2;
end

if (nargin < 7) || isempty(chopOffHighFreq),
  chopOffHighFreq=0; 
  % chop off high frequencies in order to get the same 
  % irregularity factor in the spectrum as in the data
  % not a good idea
end

if (nargin < 6) || isempty(dT),
  dT=xx(2,1)-xx(1,1);
  if ~istime
    disp(' Warning: The sampling frequency is undetermined and set to 1 Hz.')
  end
end

if (nargin < 5) || isempty(wdef)
  wdef=1; % 1=parzen window 2=hanning window
elseif ischar(wdef)
  switch lower(wdef(1))
   case 'p',
    wdef = 1; % parzen
   case 'h',
    wdef = 2; % hanning
   otherwise
    error('unknown option')
  end
end

if (nargin < 4) || isempty(plotflag)
  plotflag = 0; 
end
if (nargout == 0) && (plotflag==0)
  plotflag = 1; 
end

if (nargin<3) 
  g=[];
end

if  (isempty(g)), 
  yy=xx;
else
  yy = dat2gaus(xx,g);
end

% By using a tapered data window to smooth the data at each
% end of the record has the effect of sharpening
% the spectral window.
% NB! the resulting spectral estimate must be 
% normalized in order to correct for the loss of
% amplitude (energy) caused by the data taper.
%taper=bingham(n);
%yy(:,2)=bingham(n).*yy(:,2);
%display('****1');
%break;
yy(:,2:m)=(yy(:,2:m)-ones(n,1)*mean(yy(:,2:m)) );%/std(yy(:,2));

max_L=300;

changed_L =0;
if (nargin < 2 || isempty(L))
  L=min(n-2,ceil(4/3*max_L));
%else
  changed_L =1;
end

%if ~strcmp(method,'psd') | changed_L, 
  r=0;
  stdev=0;
  for ix=2:m
    R=dat2cov(yy(:,[1 ix]));
    r=r+R.R(:);
    stdev=stdev+R.stdev(:);
  end
  R.stdev=stdev/(m-1);
  r=r/(m-1);
  R.R=r;
  if   changed_L,
    % finding where ACF is less than 2 st. deviations.
    L=find(abs(r(1:max_L))>2*R.stdev(1:max_L))+1; % a better L value  
    L=min(L(end),max_L);
  end  
%end

if wdef==1             % Parzen window
  if changed_L,        % modify L so that hanning and Parzen give appr. the same 
                       % result
    L=min(floor(4*L/3),n-2);disp(['The default L is set to ' num2str(L) ])
  end
  v=floor(3.71*n/L);   % degrees of freedom used in chi^2 distribution
  win=parzen(2*L-1);   % Does not give negative estimates of the spectral density
  Be=2*pi*1.33/(L*dT); % bandwidth (rad/sec)
  
elseif wdef==2         % Hanning window
  if changed_L,
    disp([' The default L is set to ' num2str(L) ])
  end
  v=floor(2.67*n/L);   % degrees of freedom used in chi^2 distribution
  win=hanning(2*L-1);  % May give negative estimates of the spectral density
  Be =2*pi/(L*dT);     % bandwidth (rad/sec)
end

nf=2^nextpow2(2*L-2);  %  Interpolate the spectrum with rate = 2
nfft=2*nf;

S=createspec('freq',freqtype);
S.tr=g;
S.note=['dat2spec(',inputname(1),')'];
S.L=L;
S.S=zeros(nf+1,m-1);

% add a nugget effect to ensure that round off errors
% do not result in negative spectral estimates

r=r+nugget;
rwin=r(1:L).*win(L:(2*L-1)); 
Rper=real(fft([rwin; zeros(nfft-(2*L-1),1); rwin(L:-1:2)],nfft));

rpmin=min(Rper);
if rpmin<0
  disp(' Warning: negative spectral estimates')
end

if strcmp(freqtype,'w')
  S.w=(0:nf)'/nf*pi/dT;              % (rad/s)
  S.S=abs(Rper(1:(nf+1),1))*dT/pi;     % (m^2*s/rad) one-sided spectrum
  S.Bw=Be;
else % freqtype == f
  S.f=(0:nf)'/nf/2/dT;               % frequency Hz if dT is in seconds
  S.S=2*abs(Rper(1:(nf+1),1))*dT;      % (m^2*s) one-sided spectrum
  S.Bw=Be/(2*pi);                      % bandwidth in Hz
end


 

N=floor(nf/10);
% cutting off high frequencies 
% in this way may not be a very good idea.

if (N>3)&&chopOffHighFreq, 
  
	% The data must be Gaussian in order for this proc to be correct.
	[g0 test cmax irr]  = dat2tr(xx,'nonlinear');
  ind=(nf-4):(nf+1);
  [S,m4,m2]=specnorm(S,0);

  while (sqrt(m4/m2^2) > irr) && (ind(1)>1) 
    S.S(ind)=0;
    [S,m4,m2]=specnorm(S,0);
    ind=ind-5;
  end
  fcut=S.w(min(ind(end)+5,nf)); % cut off frequency
  if ind(1)<1,
    disp('DAT2SPEC: Error in cutting off high frequencies, try other L-values')
  end
end



%----------------------------------- 
% 
% 	Plotting the Spectral Density 
%
%-----------------------------------

if plotflag>0
 plotspec(S,plotflag)
end

