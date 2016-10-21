function [Slopes,Steep,Data]=spec2steepdat(S,Nsim,type,lev,opt,varargin)
%SPEC2STEEPDAT Simulates Lagrange waves and extracts steepness and slopes
%                
%CALL: [Slopes,Steep,Data] = spec2steepdat(S,Nsim,type,levels,options)
%
%   Slopes      = struct with fields 
%      .up      = observed slopes at the up- and 
%      .down      down-crossings of specified levels                   
%      .meanup  = average wave profiles near up- and crossings
%      .meandown    down-crossings of mean water level
%      .meanx   = at corresponding times
%      .A       = asymmetry measure by Hilbert transfor skewness
%      .lambdaAL= asymmetry measure slope ratio at mean crossing
%   Steep       = struct with fields 
%      .ffull   = full front steepness as measured by L',L'' etc
%      .bfull   = full back steepness
%      .fhalf   = half front steepness
%      .bhalf   = half back steepness
%      .lambdaN = asymmetry measure according to front/back half period 
%   Data        = struct Lagrange waves and two derivatives
%
%   S           = orbital spectrum
%   Nsim        = number of replicates in simulation (default = 1)
%   levels      = vector of standardized levels relative relative to 
%                   (default levels = [-1 0 1 2 3]*standard deviation )
%   type        = 'time', or 'space'
%   options     = struct with fields for individual replicates
%      .plotflag - 0, no plotting (default)
%                - 1, plotting of waves and cross-covariance
%                - 2, plotting of average waves
%                - 3, both the above
%
% See also: spec2ldat, ldat2lwav, wav2slope

% Tested on Matlab 8.1
% History
%   Created by GL 2008 for Skewed time waves
%   Modified by GL 2013 for new simulation program and options
%   Added front and back crest periods February 2014 by GL

if nargin<5,
    opt=simoptset;
end
if nargin>=6,  opt  = simoptset(opt,varargin{:});
end

plotflag=opt.plotflag;

if nargin<2||isempty(Nsim)
    Nsim=1;
end

if nargin<3||isempty(type)
  type='time';
end

mom=spec2mom(S);
if nargin<4||isempty(lev)
    lev=[-1 0 1 2];
end
levels=lev*sqrt(mom(1));

L0=2*pi*sqrt(mom(1)/mom(2)); %Mean upcrossing period

Nlevels = length(lev);
Slopes = struct('up',[],'down',[],...
    'meanup',[],'meandown',[],'meanx',[],...
    'levels',[],'AH',[],'lambdaAL',[]);
Slopes.date = datestr(now);
Steep = struct('ffull',[],'bfull',[],...
    'fhalf',[],'bhalf',[],'lambdaNLS',[]);
Steep.date = datestr(now);
SteepL=Steep;

for ind=1:Nlevels,
    SlopesL.up{ind}=[];
    SlopesL.down{ind}=[];
    SlopesL.up{ind}=[];
    SlopesL.down{ind}=[];
end
SlopesL.levels=levels;
dense=5;

tic
Ldata.Z=[];
Ldata.Zprim=[];
Ldata.Zbis=[];
Ldata.date=datestr(now);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate Lagrange waves with derivatives 
for k=1:10, % Generate 10 Lagrange waves and put in LData
    L=[];
    Lsmooth=[]; 
    while isempty(Lsmooth)
        [w,x]=spec2ldat(S,opt);
        [L,Lsmooth]=ldat2lwav(w,x,type,[],dense,0);
    end
    dtt=Lsmooth.t(2)-Lsmooth.t(1);
    Lprim=gradient(Lsmooth.Z,dtt);
    Lbis=gradient(Lprim,dtt);
    Ldata.t=Lsmooth.t;
    Ldata.Z=[Ldata.Z  Lsmooth.Z];
    Ldata.Zprim=[Ldata.Zprim  Lprim];    
    Ldata.Zbis=[Ldata.Zbis  Lbis];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the last example of Gaussian and Lagrange waves 

if (plotflag==1 || plotflag==3),
    ylim=2*mean(std(Wdata.Z'));
    figure(10)
    clf
    clear axis
    subplot(221)
    WIN.t=w.t;
    WIN.Z=w.Z(round(opt.Nu/2),:);
    plot(WIN.t,WIN.Z);
    axis([0 200 -1.5*ylim 1.5*ylim])
    grid on
    title('Gaussian waves')
    xlabel('t [s]')
    ylabel('Waves from model')
    subplot(222)
    plot(Lsmooth.t,Lsmooth.Z)
    title('Lagrange waves')
    xlabel('t [s]')
    grid on
    axis([0 200 -1.5*ylim 1.5*ylim])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup main simulation of Slopes and Steep
dttw=w.t(2)-w.t(1);
N0=5*floor(L0/dttw*dense); 
SlopesL.meanup=zeros(1,2*N0+1);
SlopesL.meandown=zeros(1,2*N0+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main simulations
Ndoubles=0;
LA_=0;
for n=1:Nsim,
    waitbar(n/Nsim)
    WavesLu=[];
    WavesLd=[];    
    Lsmooth=[];
    while isempty(Lsmooth), % Make sure L is generated
        Ndoubles=Ndoubles+1;
        [w,x]=spec2ldat(S,opt);
        [L,Lsmooth]=ldat2lwav(w,x,type,[],dense,0);
    end
    
    tt=Lsmooth.t;
    lmid.Z=Lsmooth.Z;
    lmid.t=tt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find skewness characteristics for original waves 
    % Find slopes at level crossings

    slopL=wav2slope(lmid,levels,dense,1);
        for ind=1:Nlevels,
            SlopesL.up{ind}=[SlopesL.up{ind}; slopL.up{ind}];
            SlopesL.down{ind}=[SlopesL.down{ind}; slopL.down{ind}];
        end
    % End finding slopes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find front and back periods
    steepL = lwav2frontback(lmid);
    SteepL.ffull=[SteepL.ffull ; steepL.ffull];
    SteepL.bfull=[SteepL.bfull ; steepL.bfull];
    SteepL.fhalf=[SteepL.fhalf ; steepL.fhalf];
    SteepL.bhalf=[SteepL.bhalf ; steepL.bhalf];
    % End finding periods
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find shape around crossings  
    Lu=lmid.Z; 
    Ld=fliplr(Lu);
    lmean=0; 
    [Lindu_0,Lncu0]=dat2crossind(Lu,lmean,'u',true);
    [Lindd_0,Lncd0]=dat2crossind(Ld,lmean,'u',true);

    if Lncu0>0,
        for k=1:Lncu0,
            if Lindu_0(k)>N0+1 && Lindu_0(k)<length(Lu)-N0+1,
                WavesLu=[WavesLu; Lu(Lindu_0(k)-N0:Lindu_0(k)+N0)];
            end
        end
    end
    if Lncd0>0,
        for k=1:Lncd0,
            if Lindd_0(k)>N0+1 && Lindd_0(k)<length(Ld)-N0+1,
                WavesLd=[WavesLd; Ld(Lindd_0(k)-N0:Lindd_0(k)+N0)];
            end
        end 
    end
    % End finding shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(WavesLu)||isempty(WavesLd)
        n=n-1;
        disp('Empty Waves')
    else % Summarize
        LHilbert=imag(hilbert(Lu));
        LA_=LA_+skew(LHilbert);

        SlopesL.meanup=((n-1)*SlopesL.meanup+mean(WavesLu))/n;
        SlopesL.meandown=((n-1)*SlopesL.meandown+mean(WavesLd))/n;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


toc
disp(['Ndoubles =' ,num2str(Ndoubles-Nsim)])

SlopesL.AH=LA_/Nsim;
SlopesL.meanx=(-N0:N0)*dttw/dense;
SlopesL.lambdaAL=(SlopesL.meandown(N0+1)-SlopesL.meandown(N0-1));
SlopesL.lambdaAL=SlopesL.lambdaAL/(SlopesL.meanup(N0+1)-SlopesL.meanup(N0-1));

Slopes.L=SlopesL;
Slopes.L.mean=lmean;
Steep=SteepL;
Steep.lambdaNLS=mean(SteepL.ffull)/mean(SteepL.bfull);
Data=Ldata;
Data.opt=opt;
Slopes=SlopesL;
Slopes.lalpha=opt.lalpha;

if plotflag==2 || plotflag==3,
    figure(20)
    clf
    clear axis
%    subplot(221)
%    plot(SlopesW.meanx,SlopesW.meanup)
%    hold on
%    grid
%    plot(SlopesW.meanx,SlopesW.meandown,'-.')
%    axis([-10 10 -ylim ylim])
%    plot([-10 10],[wmean wmean])
%    title('Average Gaussian crossings')
%
%    subplot(222)
    plot(SlopesL.meanx,SlopesL.meanup)
    hold on
    grid
    plot(SlopesL.meanx,SlopesL.meandown,'-.')
    axis([-10 10 -ylim ylim])
    plot([-10 10],[lmean lmean])
    title('Average Lagrange crossings')
end
