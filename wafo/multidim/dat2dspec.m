function [Sd,D,Sw,Fcof,Gwt,Sxy,Sxy1] = dat2dspec(xn,pos,h,nfft,nt,method,varargin)
%DAT2DSPEC Estimates the directional wave spectrum from timeseries 
% 
%  CALL: [S, D, Sw,Fcof] = dat2dspec(W,pos,h,Nfft,Nt,method,options);
%   
%         S = a spectral density structure containing:
%             S      = 2D array of the estimated directional spectrum, size Nt x Nf. 
%             w      = frequency vector  0..2*pi*Nyquistfrequency of length Nf
%             theta  = angle vector -pi..pi of length Nt 
%                      (theta = 0 -> + x-axis, theta = pi/2 -> + y-axis) 
%             h      = water depth 
%                      ( S.S=D.S.*Sw.S(:,ones(:,Nt))' )
%        D  = estimate of spreading function as function of theta and w
%        Sw = angular frequency spectrum
%      Fcof = Fourier coefficients of the spreading function D(w,theta)
%             defined as int D(w,theta)*exp(i*n*theta) dtheta. n=0:nharm
%             size nharm+1 x Nf
%
%        W  = [t w1 w2 ... wM]  a M+1 column data matrix with sampled times 
%             and values in column 1 and column 2 to M+1, respectively.
%       pos = [x y z def bfs], matrix characterizing instruments and positions, size M x 5
%             x(j), y(j),z(j) = the coordinate positon of the j-th sensor (j=1:M)
%             def(j) = identifying parameter (1,2,... or 18) for the j-th time series  
%                      (See tran for options)
%             bfs(j) = 1, if j-th time series is to be used in final estimate of 
%                         "Best" Frequency Spectra.
%                      0, if j-th is to be excluded. (At least one must be set to 1)
%        h  = water depth (default is deep water,i.e., h = inf)
%      Nfft = length of FFT  (determines the smoothing and number of
%             frequencies, Nf = Nfft/2+1) (default 256)
%        Nt = number of angles (default 101)
%    method = 'BDM'  Bayesian Directional Spectrum Estimation Method 
%             'MLM'  Maximum Likelihood Method (default)
%             'IMLM' Iterative  Maximum Likelihood Method 
%             'MEM'  Maximum Entropy Method   (slow)
%             'EMEM' Extended Maximum Entropy Method
%   options = optional options, see specoptset for details.
%
%
% Example:
%  S  = jonswap;
%  D  = spreading(linspace(-pi,pi,51),'cos2s');
%  Sd = mkdspec(S,D,1);
%  Nx = 3; Ny = 2; Nt = 2^14; dx = 10; dy = 10;dt = 0.5;
%  F  = seasim(Sd,Nx,Ny,Nt,dx,dy,dt,1,0);
%  Z  = permute(F.Z,[3 1 2]);
%  [X,Y] = meshgrid(F.x,F.y);
%  N = Nx*Ny;
%  types = repmat(sensortypeid('n'),N,1);
%  bfs   = ones(N,1);
%  pos   = [X(:),Y(:),zeros(N,1)];
%  h = inf;
%  Se = dat2dspec([F.t Z(:,:)],[pos types,bfs],h,256,101); % seasim is possibly wrong
% 
% See also  specoptset, welch_psd, welch_cpsd, tran, w2k, plotspec

%             freq = frequency vector 0..Nyquistfrequency

% References:
% Isobe M, Kondo K, Horikawa K. (1984)
% "Extension of MLM for estimating directional wave spectrum",
% In Symposium on Description and modelling of Directional Seas,
% Copenhagen, pp 1-15
%
% Young I.R. (1994)
% "On the measurement of directional spectra",
% Applied Ocean Research, Vol 16, pp 283-294
%
% Hashimoto,N. (1997) 
% "Analysis of the directional wave spectra from field data.", Advances in 
% Coastal and Ocean Engineering, Vol.3., pp.103-143
%
% Massel S.R. and Brinkman R. M. (1998)
% "On the determination od directional wave spectra 
%   for practical applications",
%  Applied Ocean Research, Vol 20, pp 357-374
%
% Michel K. Ochi (1998),
% "OCEAN WAVES, The stochastic approach",
%  OCEAN TECHNOLOGY series 6, Cambridge, chapter 7.


% Tested on: Matlab 5.3, 5.2 
% History:
% revised Jan 2008
% -Now independent of the signal toolbox
% revised pab dec2003  
%  changed name to dat2dspec
% revised pab 08.11.2002
% -changed name to dat2dspec2 in order to incorprate the specoptset option.
% revised pab 24.10.2002
% -made it more robust by estimating Hw from data if possible.
% -added method EMEM and revised the MEM code
% revised pab 09.05.2001 
%   Due to Vengatesan Venugopal plotting is now handled correctly.
% revised pab 08.05.2001
% - added optimset to fsolve for Matlab v 5.3 and higher.
% revised pab 21.06.2000
%   - added 'MEM'
%   - added 'IMLM' and checked with testsurf and testbuoy
%   - replaced inv with pinv
%   - changed order of input arguments
%   - added dflag, ftype to input
% revised pab 06.01.2000 
%   - tested  narrowbanded cases generated with testsurf and testbuoy -> OK
%   - changed size of S.S from Nf x Nt   to     Nt x Nf
%   - added the possibility that W contains timeseries of other 
%     quantities than a surface elevation
%   - added par
%   - extended pos with def and bfs 
% revised pab 28.08.99
%  updated spectral structure
% By Per A. Brodtkorb 27.03.99

% TODO % Add option 'fem'. It is not implemented
% TODO % measurements should be scaled      
% TODO % all options in specoptset should be implemented 

  
opt = specoptset('plotflag','off','dflag','mean','ftype','w',....
     'thtype','r','nharm',2,'delay',0,'gravity',[],'wdensity',[],...
     'bet',[],'igam',[],'x_axisdir',[],'y_axisdir',[],...
     'maxIter',25,'coefAbsTol',0.01,'errorTol',[],...
     'minModelOrder',1,'relax',[]);

% If just 'defaults' passed in, return the default options in Sd
if nargin==1 && nargout <= 1 && isequal(xn,'defaults')
  Sd = opt; 
  return
end
error(nargchk(2,inf,nargin)) 


if nargin>6,
  opt  = specoptset(opt,varargin{:}); 
end

xx = xn;
if nargin<2,
  error('Must have 2 inputs arguments')
end
[n m]= size(xx);
if n<m, 
  b  = m;
  m  = n;
  n  = b; 
  xx = xx.'; 
end

if n<4,  
  error('WAFO:DAT2DSPEC','The vector must have more than 4 elements!')
end

switch m
  case {1,2,3} ,
    error('Wrong dimension of input! Must at least have 4 columns! ')   
 otherwise,   % dimension OK!     
end
if any([m-1 5]~=size(pos)),  
  error('Wrong dimension of pos, must be of size MX5'),
end

Fs = 1/(xx(2,1)-xx(1,1));% sampling frequency

% Default values
%~~~~~~~~~~~~~~~
if nargin<3||isempty(h),     h      = inf; end                 % deep water is default
if nargin<4||isempty(nfft),  nfft   = min(256,n-mod(n,2)); end %make sure nfft<=n and even
if nargin<5||isempty(nt),    nt     = 101;end
if nargin<6||isempty(method),method = 'mlm';end

nfft = min(nfft+mod(nfft,2),n-mod(n,2));

dflag  = 'mean';   %'linear'; %detrending option 'none','mean','ma'
ftype  = 'w';    % options are f=S(f,phi),w=S(w,phi)
nharm  = 2;      % number of harmonics
g=[];rho=[];bet=[];igam=[];thx=[];thy=[]; % use default values in tran

% if nargout==0, 
%   opt.plotflag = 1; 
% else 
%   opt.plotflag = 0;
% end

noverlap = floor(nfft/2);

switch opt.dflag;
  case {'mean','ma','linear','none'}, dflag =opt.dflag;
end
switch opt.ftype;
  case {'w','f'}, ftype = opt.ftype;
end
switch opt.thtype;
  case {'r','d'}, thtype = opt.thtype;
end
switch opt.plotflag
  case {'none','off'},   plotflag = 0;
  case 'final', plotflag = 1;
  case 'iter',  plotflag = 2;
  otherwise
    plotflag = opt.plotflag;
end
if ~isempty(opt.noverlap),  noverlap = opt.noverlap;end
if ~isempty(opt.window),    window   = opt.window; else window = nfft; end
%if ~isempty(opt.message),   message  = opt.message; end
%if ~isempty(opt.delay),     delay    = opt.delay; end
if ~isempty(opt.nharm),     nharm    = opt.nharm; end
if ~isempty(opt.gravity),   g        = opt.gravity; end
if ~isempty(opt.wdensity),  rho      = opt.wdensity; end
if ~isempty(opt.bet),       bet      = opt.bet; end
if ~isempty(opt.igam),      igam     = opt.igam; end
if ~isempty(opt.x_axisdir), thx      = opt.x_axisdir;end
if ~isempty(opt.y_axisdir), thy      = opt.y_axisdir; end


if nt<30,  
  warning('WAFO:DAT2DSPEC',' # of angles are low =>  spectral density may be inaccurate')
end
if nharm>=nt, 
  nharm = nt-1;
  warning('WAFO:DAT2DSPEC','Nharm must be less than Nt')
end

bfs = find(pos(:,5)>.5); % indicator if it should be used to estimate spectra
if isempty(bfs),  error('Not all the elements of bfs can be zero'),end

%thetai=(linspace(0,2*pi,nt));
thetai   = linspace(-pi,pi,nt)';
nf       = nfft/2+1;


%--------------------------------------------
% Detrend the measurements before estimation
%--------------------------------------------
switch lower(dflag)
case 'none', % do nothing
case 'mean',   mn        = mean(xx(:,2:m));
               xx(:,2:m) = (xx(:,2:m)-mn(ones(n,1),:));% remove mean
case 'linear', xx(:,2:m) = detrend(xx(:,2:m));         % remove linear trend
case 'ma',     xx(:,2:m) = detrendma(xx(:,2:m),2*nfft);% remove moving average
               dflag     = 'mean'; % must change to a valid option for psd and csd
end



% Scale measurements so that they have equal standard deviation
% and zero mean.
% This is done in order to overcome possible calibration errors
% between the different measuring devices.

stdev = std(xx(:,2:m)); 
mn    = mean(xx(:,2:m));

normfact = ones(1,m-1);
def = unique(pos(:,4).');
for ix=def
  k = find(pos(:,4)==ix);
  if any(k)
    normfact(k) = stdev(k)/sqrt(mean(stdev(k).^2)); 
  end
end

xx(:,2:m) = (xx(:,2:m)-mn(ones(n,1),:))./normfact(ones(n,1),:);% remove mean


%----------------------------------------------------------------------
% finding the spectra, matrix of cross-spectra and transfer functions
%----------------------------------------------------------------------
Hw  = zeros(m-1,nf);     % a function of frequency only
Gwt = zeros(m-1,nt,nf);  % a function of frequency and direction combined.
%EE  = Gwt;

Sxy = zeros(m-1,m-1,nf);


if 1

  Sf  = zeros(m-1,nf).';
  kw=[];
  psdopt = struct('nfft',nfft,'Fs',Fs,'window',window,'overlap',noverlap,'dflag',dflag);
  %[Sf,fi] = welch_psd(xx(:,2:m),psdopt);
  %kw = w2k(2*pi*fi,0,h);
  for ix=1:m-1,
     [Sf(:,ix), fi] = welch_psd(xx(:,ix+1),psdopt);
    [Hw(ix,:),Gwt(ix,:,:),kw] = tran(2*pi*fi,thetai,pos(ix,1:3),pos(ix,4),h,g,rho,bet,igam,thx,thy,kw);
    Sxy(ix,ix,:)  = Sf(:,ix).';
    for iy=(ix+1):m-1,
      % Note welch_cpsd(x,y)*Fs/2 == conj(csd(x,y)) = csd(y,x)
      Sxy(iy,ix,:) = welch_cpsd(xx(:,ix+1),xx(:,iy+1),psdopt);
      Sxy(ix,iy,:) = conj(Sxy(iy,ix,:));
    end
  end
else % old call requires signal toolbox. Kept just in case 
  Sf  = zeros(m-1,nf).';
  kw=[];

  for ix=1:m-1,
    [Sf(:,ix), fi] = psd(xx(:,ix+1),nfft,Fs,window,noverlap,dflag);
    [Hw(ix,:),Gwt(ix,:,:),kw]=tran(2*pi*fi,thetai,pos(ix,1:3),pos(ix,4),h,g,rho,bet,igam,thx,thy,kw);
    Sxy(ix,ix,:)  = Sf(:,ix).';
    for iy=(ix+1):m-1,
      Sxy(ix,iy,:) = csd(xx(:,ix+1),xx(:,iy+1),nfft,Fs,window,noverlap,dflag);
%      Sxy(ix,iy,:) = cpsd(xx(:,ix+1),xx(:,iy+1),window,noverlap,nfft,Fs,dflag);
      Sxy(iy,ix,:) = conj(Sxy(ix,iy,:));
    end
  end

  Sf  = Sf*2/Fs;         % make sure Sf and Sxy have the correct units
  Sxy = Sxy*2/Fs;
end

if nargout>6,
  Sxy1 = Sxy*Fs/2;
end



Sf = Sf.';
%kw = kw(:)';

%-------------------------------------------------------------------------------
%Estimate frequency spectrum for the surface elevation from the bfs timeseries
%--------------------------------------------------------------------------------
SfBest = bfsspec(Sf,Hw,pos,bfs);
%---------------------------------------------------------------------------------
%Estimate the absolute value of the transfer function H(w) from the sensor spectra
%---------------------------------------------------------------------------------
Hw = hwestimate(Sf,SfBest,Hw,pos);


%------------------------------------------
%  Normalizing the matrix of cross-spectra
%-----------------------------------------

k = find(SfBest); % avoid division by zero
Sxyn = Sxy;
for ix=1:m-1
 Sxyn(ix,ix,k) = squeeze(Sxy(ix,ix,k)).'./(SfBest(k).*Hw(ix,k).^2); %1
 %Sxyn(ix,ix,k) = squeeze(Sxy(ix,ix,k)).'./(Hw(ix,k).^2); %2
 %Sxyn(ix,ix,k) = 1; %3
 %Sxx =  squeeze(Sxy(ix,ix,k)).'; %3
  for iy=(ix+1):m-1,
    Sxyn(ix,iy,k) = squeeze(Sxy(ix,iy,k)).'./(SfBest(k).*Hw(ix,k).*Hw(iy,k)); % 1
    %Sxyn(ix,iy,k) = squeeze(Sxy(ix,iy,k)).'./(Hw(ix,k).*Hw(iy,k)); % 2
    %Syy =  squeeze(Sxy(iy,iy,k)).'; %3
    %Sxyn(ix,iy,k) = squeeze(Sxy(ix,iy,k)).'./(sqrt(Sxx.*Syy)); % 3
    Sxyn(iy,ix,:) = conj(Sxyn(ix,iy,:));
  end  
end

switch lower(method)
  case {'mlm'} % maximum likelihood method
   DS = mlm(Sxyn,Gwt,thetai,fi,k,opt);
 case {'imlm'} % iterative maximum likelihood method
   DS = imlm(Sxyn,Gwt,thetai,fi,k,opt);
  case 'mem' , %Maximum entropy method
    DS = mem(Sxyn,Gwt,thetai,fi,k,opt);
  case 'emem' , %Extended Maximum entropy method      
    DS = emem(Sxyn,Gwt,thetai,fi,k,opt);
  case 'bdm', %   
    DS = bdm(Sxyn,Gwt,thetai,fi,k,opt);  
  case 'fem', % Fourier Expansion method

  error('FEM method not implemented yet')
otherwise , error('Unknown method')
end

if 0,%used for debugging
  Sxy2 = getcrossspectra(thetai,Gwt,DS);
  for ix = 1:m-1,
    for iy = ix:m-1,
      %subplot(2,1,1), plot(fi(k),squeeze(real(Sxy2(ix,iy,k))),'r',fi(k),squeeze(real(Sxyn(ix,iy,k))),'g')
      %subplot(2,1,2), plot(fi(k),squeeze(real(Sxy2(ix,iy,k))),'r',fi(k),squeeze(real(Sxyn(ix,iy,k))),'g')
      subplot(2,1,1), semilogy(fi(k),eps+abs(squeeze(real(Sxy2(ix,iy,k)))-squeeze(real(Sxyn(ix,iy,k)))),'g')
      subplot(2,1,2), semilogy(fi(k),eps+abs(squeeze(real(Sxy2(ix,iy,k)))-squeeze(real(Sxyn(ix,iy,k)))),'g')
      disp('hit any key'),pause, disp('hit any key')
    end
  end
end

Sd       = createspec('dir',ftype);
Sd.theta = thetai(:);
Sd.h     = h;
Sd.norm  = 0;
Sd.note  = ['dat2dspec(' inputname(1) ')'];

if strcmp(ftype,'w'),
  fi   = fi*2*pi;
  SfBest   = SfBest/2/pi;
  Sd.w = fi(:);
else
  Sd.f = fi(:);
end
Sd.S = (SfBest(ones(nt,1),:).*DS);
Sd = ttspec(Sd,opt.ftype,thtype);

if nargout>1
 D      = Sd;
 D.S    = DS;
 D.note = [D.note ', D(theta,w)'];

end

if nargout>2
 Sw      = createspec('freq',ftype);
 Sw.S    = SfBest(:);
 Sw.note = [D.note ', S(w)'];
 Sw.(ftype) = fi(:); 
end

if nargout>3
  %   int D(w,theta)*exp(i*n*theta) dtheta.
  if 0,
    Fcof1 = 2*pi*ifft(DS,[],1); % integration by FFT
    Pcor = [1; exp(sqrt(-1)*(1:nharm).'*thetai(1))]; % correction term to get
                                                     % the correct integration limits
    Fcof = Fcof1(1:1+nharm,:).*Pcor(:,ones(1,nf));
  else
    [a,b]=fourier(thetai,DS,2*pi,nharm);
    Fcof = (a+sqrt(-1)*b)*pi;
  end
end
%Sd = ttspec(Sd,options.ftype);
%----------------------------------- 
% 
% 	Plotting the Spectral Density 
%
%-----------------------------------
if plotflag>0,
 plotspec(Sd,plotflag) 
end
return  















