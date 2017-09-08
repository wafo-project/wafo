function W = spec2wave(Spec,options,varargin)
%SPEC2WAVE Spectral simulation of space-time Gaussian wave
%           
%CALL: W = spec2wave(spec,options);
% 
%   W    = Gaussian wave structure W with fields 
%     .Z = matrix of size [Nx Nt] 
%     .x = space coordinates, length Nxalong x-axis 
%     .t = time coordinates, length Nt
%
%   Spec = S, spectral density structure in 
%             angular frequency ('w') or frequency ('f') form 
%   options = struct with fields 
%       .Nt  = giving  Nt  time points.  (default length(S)-1=n-1).
%                If Nt>n-1 it is assummed that  S.S(k)=0  for  k>n-1
%       .Nx  = giving  Nx  space points (defult = Nt)
%       .dt  = step in grid (default dt is defined by the Nyquist freq) 
%       .dx  = step in grid (default dx is defined by the Nyquist freq)
%      (.u     = if non-empty and = [u1 u2 Nu] the u-vector will be set to 
%                u = linspace(u1,u2,Nu), ONLY TESTED for ffttype='ffttime'
%                if empty, then u = linspace(0,(Nu-1)*du,Nu))
%       .iseed  - method or starting seed for the random number generator 
%                (default = 'shuffle')
%       .ffttype - 'fftspace', fft over space, loop over time
%                   generate space series with evolvement 
%                   over time (useful if Nu > Nt),  
%                - 'ffttime', fft over time, loop over space (default) 
%                   generate time series with evolvement 
%                   over space (useful if Nt > Nu),  
%                - 'ffttwodim', 2D-fft over time and space.  
%
% Example of spec2wave 
%
%    S=jonswap; opt=simoptset; 
%    opt=simoptset(opt,'dt',0.25,'dx',0.25);
%    w = spec2wave(S,opt)
%    subplot(211)

% See also: spec2sdat, spec2field, cov2sdat, gaus2dat

% Modified 2015-05-22 by GL to allow more flexible u-vector
% Modified 2015-02-02 by GL

if nargin<2,
    opt=simoptset;
else
    opt=options;
end
if nargin>=3, opt = simoptset(opt,varargin{:});
end

if isfield(opt,'iseed') 
    iseed=opt.iseed;
else 
    iseed=[];
end

if verLessThan('matlab','7.12'),
    if isempty(iseed) || strcmp(iseed,'shuffle'),
        rand('seed',double(int32(sum(100*clock))));
    elseif isnumeric(iseed)
        rand('seed',int32(iseed));
    else
        rand('seed',double(int32(sum(100*clock))));
    end
else
    if isempty(iseed)
        iseed='shuffle';
    end
    rng(iseed); 
end

Nu=opt.Nu;
Nt=opt.Nt;
du=opt.du;
dt=opt.dt;

if strcmpi(opt.ffttype,'ffttime') || isempty(opt.ffttype')
    fftt=1;
elseif strcmpi(opt.ffttype,'fftspace')
    fftt=2;
else fftt=3;
end

S=Spec;
if ~isfield(S,'g')
    S.g=gravity;
end
if ~isfield(S,'h')
    S.h=Inf;
end
if strcmpi(S.type,'k1d')
    S=spec2spec(S,'freq'); % Make  S  to type  'freq'
end
if isfield(S,'f') % Make sure spectrum is in angular frequency
    S.w=S.f*2*pi;
    S.S=S.S/2/pi;
    S=rmfield(S,'f');
end

%U=opt.u;
%if ~isempty(U),
%    Nu=U(3);
%    u1=U(1);
%    u2=U(2);
%    u=linspace(u1,u2,Nu);
%    u=u';
%    du=(u2-u1)/Nu;
%end

if isempty(dt),
    dt=pi/S.w(end);  % Default set  dt  by the Nyqvist frequency
end
if isempty(Nt)
    Nt=2^nextpow2(length(S.w));
end
if isempty(du)
    du=pi/w2k(S.w(end),[],S.h);
end 
if isempty(Nu)
    Nu=2^nextpow2(length(S.w));
end

%***********************
% fftt = 1: fft over t and loop over u
% fftt = 2: fft over u and loop over t
% fftt = 3: 2D-fft over t and u
%***********************

if fftt==1, %timefft
    if S.w(end) > pi/dt,
        disp('Too large time step. Aliasing may occur')
        disp(['dt changed to ' num2str(pi/S.w(end))]) 
        dt=pi/S.w(end);
    end    

    S1 = specinterp(S,dt); 
    % Interpolates and pads spectrum with zeros up to 
    % correct max frequency pi/dt. It also increases 
    % the number of frequencies up to 2^m+1, including w=0. 

    S   = S1; clear S1
    if mod(Nt,2),
        Nt=Nt+1; % make sure it is even
    end 
    domega = 2*pi/(Nt*dt);
    omega  = (1:Nt)'*domega;
    kappa = w2k(omega,[],S.h);
    
elseif fftt==2, % spacefft, change spectrum to wave number type 'k'
    S = spec2spec(S,'k1d'); % Only for fft over u
    if S.k(end) > pi/du,
        disp('Too large space step. Aliasing may occur')
        disp(['du changed to ' num2str(pi/S.k(end))]) 
        du=pi/S.k(end);
    end
    
    S1 = specinterp(S,du); 
    % Interpolates and pads spectrum with zeros up to 
    % correct max wave number pi/du. It also increases 
    % the number of wavenumbers up to 2^m+1, including k=0. 

    S   = S1; clear S1
    if mod(Nu,2),
        Nu=Nu+1; % make sure it is even
    end 
    dkappa = 2*pi/(Nu*du);
    kappa = (1:Nu)'*dkappa;
    omega = k2w(kappa,[],S.h);

elseif fftt==3,
    if S.w(end) > pi/dt,    
        dt=pi/S.w(end);
        disp('Too large time step. Aliasing may occur')
        disp(['dt changed to ' num2str(dt)]) 
    end    
    if mod(Nt,2),
        Nt=Nt+1; % make sure it is even
    end 

    if w2k(S.w(end),[],S.h) > pi/du,
        du=pi/w2k(S.w(end),[],S.h);
        disp('Too large space step. Aliasing may occur')
        disp(['du changed to ' num2str(du)]) 
    end
    if mod(Nu,2),
        Nu=Nu+1; % make sure it is even
    end 
    
    S1=specinterp(S,dt);
    S=S1; clear S1;
        
    domega = 2*pi/(Nt*dt);
    omega  = (1:Nt)*domega;
    dkappa = 2*pi/(Nu*du);
    kappa = (1:Nu)'*dkappa;
end

% Prepare output variables
% ------------------------
t=(0:Nt-1)*dt; u=(0:Nu-1)*du;
w.Z=zeros(Nu,Nt);
w.x=u;
w.t=t;
w.note=Spec.note;
x.note1='Horizontal Lagrange component';

if fftt==1, %timefft
    % Interpolate spectrum to get it at correct frequencies
    % -----------------------------------------------------------
    Sinterp = interp1(S.w,S.S,omega,'pchip');
    Si=[0;Sinterp(2:end)]*domega;
    S0=spec2mom(Spec);
    testsum=sum(Si);
    if testsum<0.99*S0,
      disp('WARNING: too small dt or Nt. Spectrum information lost')
      disp(['sigma^2 = ' num2str(S0) '  Approx sigma^2 = ' num2str(testsum)]) 
    end  
    Si=max(0,Si);
    Sis=sqrt(Si); 

    % Generate standard normal random numbers for the simulations ...
    % -----------------------------------------------------------
    Zt = Sis'.*(rndnorm(0,1,1,Nt) + 1i*rndnorm(0,1,1,Nt));

    % and propagate in space, wave moving to the right (increasing u)
    % ----------------------

    Ku=u'*kappa';
    Zut = exp(1i*Ku).*repmat(Zt,Nu,1);

    % Generate w and x time series
    % ----------------------------
    w.Z=real(fft(Zut')');
    w.Z=flipud(w.Z);
   
elseif fftt==2, % spacefft
    % Interpolate spectrum to get it at correct wavenumbers
    % -----------------------------------------------------------
    Sinterp = interp1(S.k,S.S,kappa,'pchip');
    Si=[0;Sinterp(2:end)]*dkappa;
    S0=spec2mom(Spec);
    testsum=sum(Si);
    if testsum<0.99*S0,
      disp('WARNING: too small du or Nu. Spectrum information lost')
      disp(['sigma^2 = ' num2str(S0) '  Approx sigma^2 = ' num2str(testsum)]) 
    end  
    Si=max(0,Si);
    Sis=sqrt(Si); 

    % Generate standard normal random numbers for the simulations ...
    % -----------------------------------------------------------
    Zu = Sis.*(rndnorm(0,1,Nu,1) + 1i*rndnorm(0,1,Nu,1));

    % and propagate in time, wave moving to the right (increasing u)
    % ----------------------
    Wt=omega*t;
    Zut = exp(-1i*Wt).*repmat(Zu,1,Nt);
    
    % Generate w and x time/space series
    % ----------------------------
    w.Z=real(fft(Zut));
    w.Z=flipud(w.Z);
    
elseif fftt==3,
    % Interpolate spectrum to get it at correct frequencies
    % -----------------------------------------------------------
    Sinterp = interp1(S.w,S.S,omega,'pchip');
    Si=[0 Sinterp(2:end)]*domega;
    S0=spec2mom(Spec);
    testsum=sum(Si);
    if testsum<0.99*S0,
      disp('WARNING: too small dt or Nt. Spectrum information lost')
      disp(['sigma^2 = ' num2str(S0) '  Approx sigma^2 = ' num2str(testsum)]) 
    end  
    Si=max(0,Si);
    Sis=sqrt(Si); 
    
    % Generate standard normal random numbers for the simulations ...
    % -----------------------------------------------------------
    Zt = Sis.*(rndnorm(0,1,1,Nt) + 1i*rndnorm(0,1,1,Nt));
    Zut=zeros(Nu,Nt);
    for j=1:Nt,
        radindex=ceil(w2k(omega(j),[],S.h)/dkappa);
        if radindex<=Nu,
            Zut(radindex,j)=Zt(j);
        end
    end

    w.Z=real(fft2(Zut));
    w.Z=flipud(w.Z);

end 

mom=spec2mom(Spec);
w.std=std(w.Z(:));
meanper=2*pi*sqrt(mom(1)/mom(2));
w.meanperiod=meanper;
w.meanwavelength=gravity/(2*pi)*meanper^2;
W = w;



