function Slopes=spec2slopedat(S,Nsim,type,lev,options,varargin)
%SPEC2SLOPEDAT Simulates Lagrange waves and extracts slopes at crossings 
%
%CALL: Slopes = spec2slopedat(S,Nsim,type,levels,options)
%
%   Slopes   = struct array with observed slopes at the up- and 
%              down-crossings of specified levels
%
%   S        = orbital spectrum
%   Nsim     = number of replicates in simulation (default = 1)
%   levels   = vector of standardized levels relative to zero 
%                   (default = [-1 0 1 2]*standard deviation)
%   type     = 'space' (default), 'time', or 'approxtime'
%   options  = struct with fields for individual replicates
%      .Nt       = giving  Nt  time points.  (default length(S)-1=n-1).
%                  If Nt>n-1 it is assummed that S.S(k)=0 for all k>n-1
%      .Nu       = giving  Nu  space points (defult = Nt)
%      .du       = step in grid (default dt is defined by the Nyquist freq)
%      .dt       = step in grid (default dt is defined by the Nyquist freq) 
%      .lalpha   = alpha value for modified Lagrange
%      .lbeta    = beta value for modified Lagrange
%      .ffttype  = 'ffttime' (default), 'fftspace', 'ffttwodim'
%      .iseed    - setting for random number generator,  
%                  default = 'shuffle', [ int32 ]
%      .plotflag - 'off', no plotting (default)
%                - 'on' 
%
% Example:
%   S=jonswap; opt=simoptset; mom=spec2mom(S);
%   opt=simoptset(opt,'dt',0.25,'du',0.25)
%   Nsim=100;
%   levels=[0 1 2];
%   Slopes=spec2slopedat(S,Nsim,'time',levels,opt)
%
% Used by spec2lasym
% See also: spec2ldat, ldat2lwav, wav2slope

% Tested on Matlab 8.1, 8.6
% History
%   Created by GL 2008 for Skewed time waves
%   Modified by GL 2013 for new simulation program and options
%   Revised Nov 1, 2014, by GL to Matlab 8 standard
%   Modified March 2015 to use ldat2lslope for both time and space

if nargin<5,
    opt=simoptset;
elseif nargin>=6,  
    opt = simoptset(options,varargin{:});
else
    opt=options;
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

if logical(opt.plotflag>0) || opt.plotflag==true,
    plotting=true;
else
    plotting=false;
end

if nargin<2 || isempty(Nsim)
    Nsim=1;
end

if nargin<3 || isempty(type)
  type='space';
end

mom=spec2mom(S);
if nargin<4 || isempty(lev) 
    levels=[-1 0 1 2]*sqrt(mom(1));
else
    levels=lev*sqrt(mom(1));
end

[w0,x0]=spec2ldat(S,opt);
L0ti=w0.meanperiod; 
L0sp=w0.meanwavelength; 

Nindex=length(levels);
Slopes = struct('up',[],'down',[],...
    'meanwaveup',[],'meanwavedown',[],'meanwavex',[],...
    'levels',[],'A',[],'lambdaAL',[],'lalpha',opt.lalpha);
for ind=1:Nindex,
    Slopes.up{ind}=[];
    Slopes.down{ind}=[];
end
A_=0;
tic
if strcmp(type,'space'),
    tic
    dense=5;
    mwl=0;
    N0scale=max(100*floor(L0sp/50),50);
    delta=(w0.u(2)-w0.u(1))/dense;
    N0=N0scale/delta;

    mvup=[];
    mvdo=[];
    
    en=0;
    for n=1:Nsim, 
        waitbar(n/Nsim)
        
        Wavesu=[];
        Wavesd=[];
        Lsmooth=[];
        
        [w,x]=spec2ldat(S,opt);
        [~,Lsmooth]=ldat2lwav(w,x,type,[],dense);
        slop=ldat2lslope(w,x,type,levels);
        for ind=1:Nindex,
            Slopes.up{ind}=[Slopes.up{ind}; slop.up{ind}];
            Slopes.down{ind}=[Slopes.down{ind}; slop.down{ind}];
        end

        if ~isempty(Lsmooth)
            en=en+1;
            Slopes.meanwavex=(-N0:N0)*(Lsmooth.u(2)-Lsmooth.u(1)); 
            Lu=Lsmooth.Z'; 
            Ld=fliplr(Lu);
            [indu_0,ncu0]=dat2crossind(Lu,mwl,'u',true);
            [indd_0,ncd0]=dat2crossind(Ld,mwl,'u',true);

            if ncu0>0,
                for k=1:ncu0,
                    if logical(indu_0(k)>N0+1) && logical(indu_0(k)<length(Lu)-N0+1),
                        Wavesu=[Wavesu; Lu(indu_0(k)-N0:indu_0(k)+N0)];
                    end
                end 
            end
            if ncd0>0,
                for k=1:ncd0,
                    if logical(indd_0(k)>N0+1) && logical(indd_0(k)<length(Ld)-N0+1),
                        Wavesd=[Wavesd; Ld(indd_0(k)-N0:indd_0(k)+N0)];
                    end
                end 
            end
            if isempty(Wavesu) || isempty(Wavesd),
%               do nothing
            else
                LH=imag(hilbert(Lu));
                A_=A_+skewness(LH);
                mvup=[mvup; Wavesu];
                mvdo=[mvdo; Wavesd];
            end
        end
    end
    toc
    Slopes.meanwaveup=mean(mvup);
    Slopes.meanwavedown=mean(mvdo);
    disp([num2str(en) ' non-empty waves generated'])
else
    if strcmp(type,'time'),
        tic
        dense=10;
        mwl=0;
        delta=opt.dt/dense;
        N0=round(1.5*L0ti/delta); 
        mvup=[];
        mvdo=[];
        en=0;
        for n=1:Nsim,
            waitbar(n/Nsim)
            Wavesu=[];
            Wavesd=[];
            Lsmooth=[];
        
            [w,x]=spec2ldat(S,opt);
            [~,Lsmooth]=ldat2lwav(w,x,type,[],dense);
    
            slop=ldat2lslope(w,x,'time',levels);
            for ind=1:Nindex,
                Slopes.up{ind}=[Slopes.up{ind}; slop.up{ind}'];
                Slopes.down{ind}=[Slopes.down{ind}; slop.down{ind}'];
            end
            if ~isempty(Lsmooth),
                Lu=Lsmooth.Z; 
                Ld=fliplr(Lu);
                [indu_0,ncu0]=dat2crossind(Lu,mwl,'u',true);
                [indd_0,ncd0]=dat2crossind(Ld,mwl,'u',true);
        
                if ncu0>0,
                    for k=1:ncu0,
                        if indu_0(k)>N0+1 && indu_0(k)<length(Lu)-N0+1,
                            Wavesu=[Wavesu; Lu(indu_0(k)-N0:indu_0(k)+N0)];
                        end
                    end
                end

                if ncd0>0,
                    for k=1:ncd0,
                        if indd_0(k)>N0+1 && indd_0(k)<length(Ld)-N0+1,
                            Wavesd=[Wavesd; Ld(indd_0(k)-N0:indd_0(k)+N0)];
                        end
                    end 
                end
                
                if isempty(Wavesu) || isempty(Wavesd),
%                    do nothing
                else
                en=en+1;
                    LH=imag(hilbert(Lu));
                    A_=A_+skew(LH);
                    mvup=[mvup; Wavesu];
                    mvdo=[mvdo; Wavesd];
                end
            end
        end
    end
    toc
    Slopes.meanwaveup=mean(mvup);
    Slopes.meanwavedown=mean(mvdo);
    Slopes.meanwavex=(-N0:N0)*delta;%(Lsmooth.t(2)-Lsmooth.t(1)); %delta/dense;
    disp([num2str(en) ' non-empty waves generated'])
end  

Slopes.A=A_/en;
Slopes.lambdaAL=(Slopes.meanwaveup(N0+1)-Slopes.meanwaveup(N0-1));
Slopes.lambdaAL=Slopes.lambdaAL/(Slopes.meanwavedown(N0+1)-Slopes.meanwavedown(N0-1));
Slopes.note=['Slopes in ' type ' waves'];
Slopes.levels=levels;
