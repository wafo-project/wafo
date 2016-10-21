function [TTfu,varargout]=spec2timeslopepdf(S,y,levels,opt,varargin) 
%SPEC2TIMESLOPEPDF Computes pdf for slopes at crossings of time waves 
%
%CALL: [TTfu,TTfd,STfu,STfd,VTfu,VTfd] = ... 
%                       spec2timeslopepdf(S,y,levels,options,varargin)
%
%   XTfu, XTfd = pdf for slopes at up- and down-crossings of levels 
%                according to Sections 5.1.1, 5.1.2, 5.1.3 in 
%                Adv Appl Probab 42 (2010) 489-508. 
%                TT = time slopes at time crossings
%                ST = space slopes at time crossings
%                VT = velocity at time crossings
%                Any number of output cdf:s can be specified, 
%                starting with TTFu, [with optional TTFd, ...]      
%                pdf is computed from a smoothed gradient of the
%                simulated cdf; see spec2timeslopecdf 
%
%                For space waves, use spec2spaceslopepdf 
%
%   y         = pdf calculated at  y
%               Note: the accuracy will depend on the y-spacing
%               and on an interior smoothing parameter  smoothp  
%   levels    = vector of relaive levels 
%               (default = [-1 0 1 2])
%   options   = struct with fields (plus some more)
%       .lalpha = alpha value for modified Lagrange
%       .lbeta  = beta value for modified Lagrange
%
% Example:
%   S=jonswap; mom=spec2mom(S);
%   opt=simoptset('du',0.125,'Nu',256,'dt',0.125);
%   levels=[0 1 2]*sqrt(mom(1));
%   y=linspace(0,8,101);
%   [TTfu,TTfd]=spec2timeslopepdf(S,y,levels,opt,'lalpha',1)
%   clf
%   plot(TTfu.x,TTfu.f{1},TTfu.x,TTfu.f{2},TTfu.x,TTfu.f{3}); hold on
%   plot(TTfd.x,TTfd.f{1},'-.',TTfd.x,TTfd.f{2},'-.',TTfd.x,TTfd.f{3},'-.')
%   title('Slope PDF at up- and downcrossings, asymmetric time waves')
%   axis([0 8 0 1])
%
% See also: spec2ldat, spec2slopepdf, ldat2lwav, wav2slope

% Used in Adv Appl Probab 42 (2010) 489-508. 
% Modified November 2014 for use in WafoL

nout=max(nargout,1)-1;

smoothp=1;%0.99999;
Nsim=300; % Decrease/increase for increased speed/accuracy

m=spec2mom(S);
if nargin<4,
    opt=simoptset;
end

if nargin>=5,  opt  = simoptset(opt,varargin{:});
end

alpha=opt.lalpha;
beta=opt.lbeta;

if nargin<2 || isempty(y),
    y=(-4:0.01:4)*sqrt(m(2));
end

if nargin<3 || isempty(levels) 
    mom=spec2mom(S);
    levels=[-1 0 1 2 3]*sqrt(mom(1));
end

Nindex=length(levels);
TTFu = struct('x',y,'levels',levels);
if nout > 0, TTFd = struct('x',y,'levels',levels); end 
if nout > 1, STFu = struct('x',y,'levels',levels); end
if nout > 2, STFd = struct('x',y,'levels',levels); end
if nout > 3, VTFu = struct('x',y,'levels',levels); end
if nout > 4, VTFd = struct('x',y,'levels',levels); end
TTfu = struct('x',y,'levels',levels);
dx=y(2)-y(1);
if nout > 0, TTfd = struct('x',y,'levels',levels); end 
if nout > 1, STfu = struct('x',y,'levels',levels); end
if nout > 2, STfd = struct('x',y,'levels',levels); end
if nout > 3, VTfu = struct('x',y,'levels',levels); end
if nout > 4, VTfd = struct('x',y,'levels',levels); end
for ind=1:Nindex,
    TTFu.f{ind}=[];
    if nout > 0, TTFd.f{ind}=[]; end
    if nout > 1, STFu.f{ind}=[]; end
    if nout > 2, STFd.f{ind}=[]; end
    if nout > 3, VTFu.f{ind}=[]; end
    if nout > 4, VTFd.f{ind}=[]; end
end
type = 'time'; % Routine serves only type='time'
if strncmp(type,'time',1);
    [rww00,rwx00,rxx00]=spec2lcov(S,0,0,1,alpha,beta);
    [rwwt0,rwxt0,rxxt0]=spec2lcov(S,0,0,2,alpha,beta);
    [rwwu0,rwxu0,rxxu0]=spec2lcov(S,0,0,3,alpha,beta);
    [rwwtt,rwxtt,rxxtt]=spec2lcov(S,0,0,4,alpha,beta);
    [rwwuu,rwxuu,rxxuu]=spec2lcov(S,0,0,5,alpha,beta);
    [rwwtu,rwxtu,rxxtu]=spec2lcov(S,0,0,6,alpha,beta);
    SZZ=[rww00.R rwx00.R; rwx00.R rxx00.R];
    SYZ=[rwwt0.R rwxt0.R; rwwu0.R rwxu0.R;...
        -rwxt0.R rxxt0.R; -rwxu0.R rxxu0.R];
    SYY=[rwwtt.R rwwtu.R rwxtt.R rwxtu.R;...
         rwwtu.R rwwuu.R rwxtu.R rwxuu.R;...
         rwxtt.R rwxtu.R rxxtt.R rxxtu.R;...
         rwxtu.R rwxuu.R rxxtu.R rxxuu.R];
    SY_Z=SYY-SYZ*(SZZ\SYZ');
    CH=chol(SY_Z);
    sX=sqrt(rxx00.R);
%    uu=sX*(-4:0.01:4)
    uu=sX*(-4:0.02:4);
    yy=y;
    yy=[yy inf];
    Nu=length(uu);
    Ny=length(yy);
    NC=125; % function g_v simulates NC 4-dimensional vectors at a time
    TTF_up=zeros(Nu,Ny);
    if nout > 0, TTF_do=zeros(Nu,Ny); end
    if nout > 1, STF_up=zeros(Nu,Ny); end
    if nout > 2, STF_do=zeros(Nu,Ny); end
    if nout > 3, VTF_up=zeros(Nu,Ny); end
    if nout > 4, VTF_do=zeros(Nu,Ny); end
    tic    
for vv=1:Nindex,
    v=levels(vv);

    for kk=1:Nu,
        u=uu(kk);
        TTIndsum_up=zeros(1,Ny);
        if nout > 0, TTIndsum_down=zeros(1,Ny); end
        if nout > 1, STIndsum_up=zeros(1,Ny); end
        if nout > 2, STIndsum_down=zeros(1,Ny); end
        if nout > 3, VTIndsum_up=zeros(1,Ny); end
        if nout > 4, VTIndsum_down=zeros(1,Ny); end
        
        for N=1:Nsim,
            [TTInd_up,TTInd_down,...
                STInd_up,STInd_down,...
                VTInd_up,VTInd_down]=g_y(v,yy,u,SYZ,SZZ,CH,NC,nout);
            TTIndsum_up=TTIndsum_up+TTInd_up;
            if nout > 0, TTIndsum_down=TTIndsum_down+TTInd_down; end
            if nout > 1, STIndsum_up=STIndsum_up+STInd_up; end 
            if nout > 2, STIndsum_down=STIndsum_down+STInd_down; end
            if nout > 3, VTIndsum_up=VTIndsum_up+VTInd_up; end 
            if nout > 4, VTIndsum_down=VTIndsum_down+VTInd_down; end 
        end  % end for N=1:Nsim,
        
        TTInd_up=TTIndsum_up/Nsim;
        TTF_up(kk,:)=TTInd_up;        
        if nout > 0, 
            TTInd_down=TTIndsum_down/Nsim; 
            TTF_do(kk,:)=TTInd_down;        
        end
        if nout > 1, 
            STInd_up=STIndsum_up/Nsim; 
            STF_up(kk,:)=STInd_up;            
        end
        if nout > 2, 
            STInd_down=STIndsum_down/Nsim; 
            STF_do(kk,:)=STInd_down;        
        end
        if nout > 3, 
            VTInd_up=VTIndsum_up/Nsim; 
            VTF_up(kk,:)=VTInd_up;        
        end
        if nout > 4, 
            VTF_do(kk,:)=VTInd_down;
            VTInd_down=VTIndsum_down/Nsim; 
        end 
        
    end % end  for kk=1:nu, 
    toc
    
%    wxpdf=mvnpdf([v 0],[zeros(Nu,1) uu'],SZZ)*ones(1,Ny);
    wxpdf=pdfnorm2d(repmat([v 0],Nu,1),[zeros(Nu,1) uu'],SZZ)*ones(1,Ny);
    TTFUP=simpson(uu',wxpdf.*TTF_up);
    TTFUP=TTFUP/TTFUP(end);   
    TTFu.f{vv}=TTFUP(1:end-1); 
    TTfu.f{vv}=max(0,cssmooth(y,gradient(TTFUP(1:end-1),dx),smoothp,y));
    
    if nout > 0, 
        TTFDO=simpson(uu',wxpdf.*TTF_do);
        TTFDO=TTFDO/TTFDO(end);  
        TTFd.f{vv}=TTFDO(1:end-1);  
        TTfd.f{vv}=max(0,cssmooth(y,gradient(TTFDO(1:end-1),dx),smoothp,y));
    end
    if nout > 1, 
        STFUP=simpson(uu',wxpdf.*STF_up);
        STFUP=STFUP/STFUP(end);
        STFu.f{vv}=STFUP(1:end-1);
        STfu.f{vv}=max(0,cssmooth(y,gradient(STFUP(1:end-1),dx),smoothp,y));
    end
    if nout > 2,
        STFDO=simpson(uu',wxpdf.*STF_do);
        STFDO=STFDO/STFDO(end);    
        STFd.f{vv}=STFDO(1:end-1);
        STfd.f{vv}=max(0,cssmooth(y,gradient(STFDO(1:end-1),dx),smoothp,y));
    end
    if nout > 3,
        VTFUP=simpson(uu',wxpdf.*VTF_up);
        VTFUP=VTFUP/VTFUP(end);
        VTFu.f{vv}=VTFUP(1:end-1);
        VTfu.f{vv}=max(0,cssmooth(y,gradient(VTFUP(1:end-1),dx),smoothp,y));
    end
    if nout > 4,
        VTFDO=simpson(uu',wxpdf.*VTF_do);
        VTFDO=VTFDO/VTFDO(end);
        VTFd.f{vv}=VTFDO(1:end-1);
        VTfu.f{vv}=max(0,cssmooth(y,gradient(VTFUP(1:end-1),dx),smoothp,y));
    end
end

if nout > 0, varargout(1)={TTfd}; end
if nout > 1, varargout(2)={STfu}; end
if nout > 2, varargout(3)={STfd}; end
if nout > 3, varargout(4)={VTfu}; end
if nout > 4, varargout(5)={VTfd}; end
end

function [TTInd_up,TTInd_down,...
          STInd_up,STInd_down, ...
          VTInd_up,VTInd_down]=g_y(v,yyy,u,SYZ,SZZ,CH,NC,nout)
    noll=[0 0 0 0];
    mu=[0 0 0 1]' + SYZ*(SZZ\[v -u]');
    mu=mu';
    Ident=eye(4); 
%    U=mvnrnd(noll,Ident,NC);
    U=rndnormnd(noll,Ident,NC);
    Y=U;
    Y=(CH'*Y)';
    Y=ones(NC,1)*mu+Y;
    TTder=Y(:,1)-Y(:,2).*Y(:,3)./Y(:,4);
    STder=Y(:,2)./Y(:,4);
    VTder=Y(:,3);
    det=abs(Y(:,1).*Y(:,4)-Y(:,2).*Y(:,3));
    Nyyy=length(yyy);
    dety=det*ones(1,Nyyy);
    
    TTdery_up=(TTder*ones(1,Nyyy))<(ones(NC,1)*yyy);
    TTdery_up=TTdery_up.*(TTder*ones(1,Nyyy)>0);
    TTInd_up=mean(dety.*TTdery_up); 
    
    if nout > 0,
        TTdery_do=(TTder*ones(1,Nyyy))>(-ones(NC,1)*yyy);
        TTdery_do=TTdery_do.*(TTder*ones(1,Nyyy)<0);
        TTInd_down=mean(dety.*TTdery_do);
    else TTInd_down=[];
    end
    
    if nout > 1, 
        STdery_up=(STder*ones(1,Nyyy))<(ones(NC,1)*yyy);
        STdery_up=STdery_up.*(TTder*ones(1,Nyyy)>0);
        STInd_up=mean(dety.*STdery_up);
    else STInd_up=[]; 
    end

    if nout > 2, 
        STdery_do=(STder*ones(1,Nyyy))<(ones(NC,1)*yyy);
        STdery_do=STdery_do.*(TTder*ones(1,Nyyy)<0);    
        STInd_down=mean(dety.*STdery_do);
    else STInd_down=[];
    end
    
    if nout > 3,
        VTdery_up=(VTder*ones(1,Nyyy))<(ones(NC,1)*yyy);
        VTdery_up=VTdery_up.*(TTder*ones(1,Nyyy)>0);
        VTInd_up=mean(dety.*VTdery_up);   
    else VTInd_up=[];
    end

    if nout > 4,
        VTdery_do=(VTder*ones(1,Nyyy))<(ones(NC,1)*yyy);
        VTdery_do=VTdery_do.*(TTder*ones(1,Nyyy)<0);
        VTInd_down=mean(dety.*VTdery_do);
    else VTInd_down=[];
    end

