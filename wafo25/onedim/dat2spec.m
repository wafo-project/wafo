function [S,fcut] = dat2spec(xn,varargin)
%DAT2SPEC Estimate one-sided spectral density from data.
% 
% CALL:  S = dat2spec(x,L,g,plotflag,p,method,dflag,ftype)
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
%             CI   = lower and upper confidence constant  
%             p    = confidence level. (Default 0.95).
%             Bw   = Bandwidth of the smoothing window which is used 
%                    in the estimated spectrum. (rad/sec or Hz)
%            
%         x =  m column data matrix with sampled times in the first column
%              and values the next columns.    
%
%         L = maximum lag size of the window function. 
%             If no value is given the lag size is set to
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
%  Method   = 'cov'   Frequency smoothing using a parzen window function
%                     on the estimated autocovariance function.  (default)
%             'psd'   Welch's averaged periodogram method with no overlapping 
%                     batches
%             'psdo'  Welch's averaged periodogram method with overlapping 
%                     batches
%             'pmem'  Maximum Entropy Method (psd using the Yule-Walker 
%                     AR method)
%             'pburg' Burg's method.
%
%  dflag    = specifies a detrending performed on the signal before estimation.
%             'mean','linear' or 'ma' (= moving average)  (default 'mean')   
%  ftype    = frequency type, 'w' or 'f'  (default 'w')
%
% Method == 'cov','psd':
%  As L decreases the estimate becomes smoother and Bw increases. If we
%  want to resolve peaks in S which is Bf (Hz or rad/sec) apart then Bw < Bf. 
% 
% Method == 'pmem','pburg':
%  L denotes the order of the AR (AutoRegressive) model. 
%
%  NOTE: The strings method,dflag and ftype may be given anywhere after x 
%        and in any order.
%
% Example
%  x = load('sea.dat');
%  S = dat2spec(x);
%  plotspec(S)
% 
% See also   dat2tr, dat2cov

% Secret option: if chopOffHighFreq==1
%  Chop off high frequencies in order to 
%  get the same irregularity factor in the spectrum as 
%  in x.


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

%

% Tested on: Matlab 5.3
% history:
% revised pab 03.11.2000
% - changed call from chi2inv to wchi2inv 
% revised pab 10.07.00
% - fixed a bug for m>2 and g is given: replaced dat2gaus call with tranproc
% revised pab 12.06.2000
% - added method,dflag,ftype to input arguments 
% - removed dT wdef from input arguments
% revised pab 14.02.2000
%  - added rate for interpolation in frequency domain
% revised pab 20.01.2000
%  - added detrending option linear, ma mean
%  - added tapering of data before estimation
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


% Initialize constants 
%~~~~~~~~~~~~~~~~~~~~~
nugget   = 0; %10^-12;
rate     = 2; % interpolationrate for frequency
tapery   = 0; % taper the data before the analysis
wdef     = 1; % 1=parzen window 2=hanning window, 3= bartlett window

% Default values:
%~~~~~~~~~~~~~~~~
L        = []; 
g        = [];
plotflag = 0;
p        = 0.95;
chopOffHighFreq=0;   % chop off high frequencies in order to get the same 
                     % irregularity factor in the spectrum as in the data
                     % may not be a good idea => default is 0
		     
method   = 'cov';  % cov. other options from signal toolbox: psd = welch's method, pyulear,pmem =
                   % maximum entropy method 
dflag    = 'mean'; %'ma','linear' 'mean','none' % detrending option
ftype    = 'w'  ;  %options are 'f' and 'w'


P  = varargin;
Np = length(P);
if Np>0
  strix = zeros(1,Np);
  for ix=1:Np, % finding symbol strings 
    strix(ix)=ischar(P{ix});
  end
  k = find(strix);
  Nk = length(k);
  if Nk>0
    if Nk>3,      
      warning('WAFO:DAT2SPEC','More than 3 strings are not allowed'),
    end
    for ix = k
      switch lower(P{ix})
        case {'f','w'},                                 ftype  = P{ix};
        case {'cov','pmem','mem','psd','psdo','pburg'}, method = P{ix};
        case {'mean','ma','linear','none'},              dflag = P{ix};
        otherwise,
          warning('WAFO:DAT2SPEC',['Unknown option:' P{ix}])
      end
    end
    Np = Np-Nk;
    P  = {P{find(~strix)}}; % remove strings from input
  end

  if Np>0 && ~isempty(P{1}), L        = P{1};end
  if Np>1 && ~isempty(P{2}), g        = P{2};end
  if Np>2 && ~isempty(P{3}), plotflag = P{3};end
  if Np>3 && ~isempty(P{4}), p        = P{4};end
  if Np>4 && ~isempty(P{5}), chopOffHighFreq = P{5};end
end  
  

%[L,g,plotflag,p,method,dflag,ftype,chopOffHighFreq] = d2schk(varargin);


if (nargout == 0) && (plotflag==0)
  plotflag = 1; 
end

xx     = xn;
[n m]  = size(xx);

if min(m,n)==1, 
  xx = [ (1:n)' xx(:)];
  n  = max(m,n);
  m  = 2;
  disp('Warning: The sampling frequency is undetermined and set to 1 Hz.')
end
dT = xx(2,1)-xx(1,1);
L = min(L,n);

if  (isempty(g)), 
  yy = xx;
else
  yy = xx;
  for ix=2:m
    yy(:,ix)= tranproc(xx(:,ix),g);
  end
  %yy = dat2gaus(xx,g);
end

%display('****1');
%break;
switch lower(dflag)
  case 'mean',
    ma        = mean(yy(:,2:m));
    yy(:,2:m) = (yy(:,2:m)-ma(ones(n,1),:) );
  case 'linear',
    yy(:,2:m) = detrend(yy(:,2:m),1); % signal toolbox detrend
  case 'ma',
    dL        = ceil(1200/2/dT); % approximately 20 min. moving average    
    yy(:,2:m) = detrendma(yy(:,2:m),dL);
    dflag     ='mean';
end
% By using a tapered data window to smooth the data at each
% end of the record has the effect of sharpening
% the spectral window.
% NB! the resulting spectral estimate must be 
% normalized in order to correct for the loss of
% amplitude (energy) caused by the data taper.
sa = std(yy(:,2:m));
if tapery
 taper     = bingham(n);
 yy(:,2:m) = taper(:,ones(1,m-1)).*yy(:,2:m);
end

max_L     = min(300,n); % maximum lag if L is undetermined
changed_L = 0;
if isempty(L)
  L = min(n-2,ceil(4/3*max_L));
  changed_L = 1;
end

if strcmp(method,'cov') || changed_L, 
  r=0;
  stdev=0;
  for ix=2:m
    R = dat2cov(yy(:,[1 ix]));
    r = r+R.R(:);
    stdev = stdev+R.stdev(:);
  end
  r       = r/(m-1);
  R.stdev = mean(sa.^2)/r(1)*stdev/(m-1);
  r       = r*mean(sa.^2)/r(1);
  R.R     = r;
  %covplot(R,150)
  if   changed_L,
    %finding where ACF is less than 2 st. deviations.
    L = find(abs(r(1:max_L))>2*R.stdev(1:max_L))+1; % a better L value  
    L = min(L(end),max_L);
    
    if wdef==1   % modify L so that hanning and Parzen give appr. the same result
      L = min(floor(4*L/3),n-2);
    end
    disp(['The default L is set to ' num2str(L) ])
  end  
end

if wdef==1             % Parzen window
  v   = floor(3.71*n/L);   % degrees of freedom used in chi^2 distribution
  win = parzen(2*L-1);     % Does not give negative estimates of the spectral density
  Be  = 2*pi*1.33/(L*dT);  % bandwidth (rad/sec)
  
elseif wdef==2         % Hanning window
  v   = floor(2.67*n/L);   % degrees of freedom used in chi^2 distribution
  win = hanning(2*L-1);    % May give negative estimates of the spectral density
  Be  = 2*pi/(L*dT);       % bandwidth (rad/sec)

else wdef==3             % Bartlett window
  v   = floor(3*n/L);   % degrees of freedom used in chi^2 distribution
     win = bartlett(2*L-1);
  Be  = 2*pi*1.33/(L*dT);  % bandwidth (rad/sec)
end

nf   = rate*2^nextpow2(2*L-2);  %  Interpolate the spectrum with rate 
nfft = 2*nf;

S      = createspec('freq',ftype);
S.tr   = g;
S.note = ['dat2spec(',inputname(1),'), Method = ' method ];
S.norm = 0; % not normalized
S.L    = L;



S.S = zeros(nf+1,m-1);

switch lower(method)
  case {'psd','psdo'} % from signal toolbox
    noverlap   = 0;
    if method(end)=='o',
      noverlap   = floor(L/2);
    end
    S.noverlap = noverlap;
    if 1,
      vararg = cell(1,3*~isempty(p));
      [Rper vararg{:}]=welch_psd(yy(:,2:m),'nfft',nfft,'window',win,'overlap',noverlap,'p',p,'dflag',dflag);
      if m>2
        Rper = mean(Rper,2);
      end
      if ~isempty(p)        
        Sc = [mean(vararg{3}(:,1:m),2),mean(vararg{3}(:,m+1:end),2)] ;
        [maxS ind] = max(Rper);
        S.CI = Sc(ind,:)/Rper(ind);
        S.p  = p;
      end
    else
      % from signal toolbox
      [Rper,Sc f] = psd(yy(:,ix),nfft,1/dT,win,noverlap,p,dflag);
      for ix=3:m
        Rper = Rper + psd(yy(:,ix),nfft,1/dT,win,noverlap,p,dflag);
      end
      Rper = Rper/(m-1);
    
     if (  ~isempty(p) ),       % Confidence interval constants
       [maxS ind] = max(Rper);
       S.CI = Sc(ind,:)/Rper(ind);
       S.p  = p;
     end
    end
  case {'pyulear','pmem','mem'} % from signal toolbox
    Rper = pyulear(yy(:,2),L,nfft,1/dT);
    for ix=3:m
      Rper = Rper + pyulear(yy(:,ix),L,nfft,1/dT);
    end
    Rper = Rper/(m-1);
  case 'pburg',  % from signal toolbox
    Rper = pburg(yy(:,2),L,nfft,1/dT);
    for ix=3:m
      Rper = Rper + pburg(yy(:,ix),L,nfft,1/dT);
    end
    Rper = Rper/(m-1);
  otherwise, % cov method
  % add a nugget effect to ensure that round off errors
  % do not result in negative spectral estimates
  r    = r+nugget;
  rwin = r(1:L).*win(L:(2*L-1)); 
  Rper = real(fft([rwin; zeros(nfft-(2*L-1),1); rwin(L:-1:2)],nfft));
  %f    =  [0:(nf)]'/nf/(2*dT);
  if (  ~isempty(p) ),
    alpha = (1-p);
    % Confidence interval constants
    S.CI = [v/invchi2( 1-alpha/2 ,v) v/invchi2( alpha/2 ,v)];
    S.p  = p;
  end
end

ind = find(Rper<0);
if any(ind)
  Rper(ind) = 0; % set negative values to zero
  warning('WAFO:DAT2SPEC','negative spectral estimates')
end

if strcmp(ftype,'w')
  S.w  = (0:nf)'/nf*pi/dT;           % (rad/s)
  S.S  = real(Rper(1:(nf+1),1))*dT/pi; % (m^2*s/rad)one sided spectrum
  S.Bw = Be;
else % ftype == f
  S.f  = (0:nf)'/nf/2/dT;            % frequency Hz if dT is in seconds
  S.S  = 2*real(Rper(1:(nf+1),1))*dT;  % (m^2*s) one sided spectrum
  S.Bw = Be/(2*pi);                    % bandwidth in Hz
end

 

N = floor(nf/10);
% cutting off high frequencies 
% in this way may not be a very good idea.

if ((N>3) && chopOffHighFreq), 
  % The data must be Gaussian in order for this proc to be correct.
  [g0 test cmax irr]  = dat2tr(xx,'nonlinear');
  ind=(nf-4):(nf+1);
  
  [Sn,m4] = specnorm(S,0);

  while (sqrt(m4) > irr) && (ind(1)>1) 
    S.S(ind)   = 0;
    [Sn,m4]    = specnorm(S,0);
    ind        = ind-5;
  end
  fcut = S.w(min(ind(end)+5,nf)); % cut off frequency
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

return

function [L,g,plotflag,p,method,dflag,ftype,chopOffHighFreq] = d2schk(P)
% D2SCHK Helper function for dat2spec.
%
% CALL  [L,g,plotflag,p,method,dflag,ftype,chopOffHighFreq]=d2schk(P) 
%
%   P = the cell array P of input arguments (between 0 and 7 elements)
%  xx = must be a two column vector.
%

% Default values:
%~~~~~~~~~~~~~~~~
L        = []; 
g        = [];
plotflag = 0;
p        = 0.95;
chopOffHighFreq=0;   % chop off high frequencies in order to get the same 
                     % irregularity factor in the spectrum as in the data
                     % may not be a good idea => default is 0
		     
method   = 'cov';  % cov. other options from signal toolbox: psd = welch's method, pyulear,pmem =
                   % maximum entropy method 
dflag    = 'mean'; %'ma','linear' 'mean','none' % detrending option
ftype    = 'w'  ;  %options are 'f' and 'w'



Np=length(P);
strix=zeros(1,Np);
for ix=1:Np, % finding symbol strings 
 strix(ix)=ischar(P{ix});
end
k = find(strix);
if any(k)
  Nk=length(k);
  if Nk>3
    warning('WAFO:DAT2SPEC','More than 3 strings are not allowed in ')
  end
  for ix = k
    switch lower(P{ix})
    case {'f','w'},                                 ftype  = P{ix};
    case {'cov','pmem','mem','psd','psdo','pburg'}, method = P{ix};
    case {'mean','ma','linear','none'},              dflag = P{ix};
    otherwise,                    
      warning('WAFO:DAT2SPEC',['Unknown option:' P{ix}])
    end
  end
  Np=Np-Nk;
  P={P{find(~strix)}}; % remove strings from input
end

if Np>0 && ~isempty(P{1}), L        = P{1};end
if Np>1 && ~isempty(P{2}), g        = P{2};end
if Np>2 && ~isempty(P{3}), plotflag = P{3};end
if Np>3 && ~isempty(P{4}), p        = P{4};end
if Np>4 && ~isempty(P{5}), chopOffHighFreq = P{5};end


return
