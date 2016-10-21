function [TTFu,varargout]=spec2timeslopecdf(S,y,lev,opt,varargin)
%SPEC2TIMESLOPECDF Computes cdf for slopes at crossings of time waves 
%
%CALL: [TTFu,TTFd,STFu,STFd,VTFu,VTFd] = ... 
%                       spec2timeslopecdf(S,y,levels,options,varargin)
%
%   XTFu, XTFd  = cdf for slopes at up- and down-crossings of levels 
%                 according to Sections 5.1.1, 5.1.2, 5.1.3 in 
%                 Adv Appl Probab 42 (2010) 489-508. 
%                 TT = time slopes at time crossings
%                 ST = space slopes at time crossings
%                 VT = velocity at time crossings
%                 Any number of output cdf:s can be specified, 
%                 starting with TTFu, [with optional TTFd, ...]   
%                 Note: Crossings with retrograd x-movement not included. 
%
%                 For space waves, use spec2spaceslopecdf 
%
%   y           = cdf calculated at  y
%   levels      = vector of relative levels (default = [-1 0 1 2])
%   options     = struct with fields (plus some more)
%       .lalpha   = alpha value for modified Lagrange
%       .lbeta    = beta value for modified Lagrange
%       .plotflag - 'off', no plotting (default)
%                 - 'on' 
%
% Example:
%   S=jonswap; mom=spec2mom(S);
%   opt=simoptset('du',0.125,'Nu',512,'dt',0.125);
%   levels=[0 1 2]*sqrt(mom(1));
%   y=linspace(0,8,1001);
%   [TTFu,TTFd]=spec2timeslopecdf(S,y,levels,opt,'lalpha',1)
%   clf
%   plot(TTFu.x,TTFu.f{1},TTFu.x,TTFu.f{2},TTFu.x,TTFu.f{3}); hold on
%   plot(TTFd.x,TTFd.f{1},'-.',TTFd.x,TTFd.f{2},'-.',TTFd.x,TTFd.f{3},'-.')
%   title('Slope CDF at up- and downcrossings, asymmetric time waves')
%   axis([0 8 0 1])
%
% See also: spec2ldat, spec2slopepdf, ldat2lwav, wav2slope

% Used in Adv Appl Probab 42 (2010) 489-508. 

% Tested in Matlab 8.1, 8.6
% History:
% Created by GL for AAP paper
% Modified November 2014 for use in WafoL

nout=max(nargout,1)-1;

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

mom=spec2mom(S);
if nargin<3 || isempty(lev) 
    levels=[-1 0 1 2 ]*sqrt(mom(1));
else
    levels=lev*sqrt(mom(1));
end

Nlevels=length(levels);
TTFu = struct('x',y,'levels',levels);
if nout > 0, TTFd = struct('x',y,'levels',levels); end 
if nout > 1, STFu = struct('x',y,'levels',levels); end
if nout > 2, STFd = struct('x',y,'levels',levels); end
if nout > 3, VTFu = struct('x',y,'levels',levels); end
if nout > 4, VTFd = struct('x',y,'levels',levels); end
for ind=1:Nlevels,
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
     SZZinv=inv(SZZ);
    SY_Z=SYY-SYZ*(SZZ\SYZ');
    CH=chol(SY_Z);
    sX=sqrt(rxx00.R);
    yy=y;
    yy=[yy inf];
%    uu=sX*(-4:0.01:4);
%    Nu=length(uu);
    Ny=length(yy);
    NC=100; % function g_v simulates NC 4-dimensional vectors at a time
    tic   
for vv=1:Nlevels,
    v=levels(vv);
    SX_W=SZZ(2,2)-SZZ(1,2)^2/SZZ(1,1);
    MX_v=SZZ(1,2)/SZZ(1,1)*v;
%    uu=sX*[-4:0.01:4]; 
    uu=MX_v+sX*(-8:0.01:8);
%    Cond_std_xu=sqrt(SY_Z(4,4));
%    umax=(invnorm(0.9999,0,1)*Cond_std_xu-1+rwxu0.R*v*SZZinv(1,2))/rwxu0.R/SZZinv(1,1);
%    uu=linspace(-4*sX,umax,500);
    Nu=length(uu);
    TTF_up=zeros(Nu,Ny);
    if nout > 0, TTF_do=zeros(Nu,Ny); end
    if nout > 1, STF_up=zeros(Nu,Ny); end
    if nout > 2, STF_do=zeros(Nu,Ny); end
    if nout > 3, VTF_up=zeros(Nu,Ny); end
    if nout > 4, VTF_do=zeros(Nu,Ny); end

    for kk=1:Nu,
        waitbar(kk/Nu);
        u=uu(kk);
        TTIndsum_upp=zeros(1,Ny);
        TTIndsum_upm=zeros(1,Ny);
        if nout > 0, 
            TTIndsum_downp=zeros(1,Ny);
            TTIndsum_downm=zeros(1,Ny);
        end
        if nout > 1, STIndsum_up=zeros(1,Ny); end
        if nout > 2, STIndsum_down=zeros(1,Ny); end
        if nout > 3, VTIndsum_up=zeros(1,Ny); end
        if nout > 4, VTIndsum_down=zeros(1,Ny); end
        Nsim=10;
        
        for N=1:Nsim,
            [TTInd_upp,TTInd_upm,TTInd_downp,TTInd_downm,...
                STInd_up,STInd_down,...
                VTInd_up,VTInd_down]=g_y(v,yy,u,SYZ,SZZ,CH,NC,nout);
            TTIndsum_upp=TTIndsum_upp+TTInd_upp;
            TTIndsum_upm=TTIndsum_upm+TTInd_upm;
            if nout > 0, 
                TTIndsum_downp=TTIndsum_downp+TTInd_downp; 
                TTIndsum_downm=TTIndsum_downm+TTInd_downm;
            end
            if nout > 1, STIndsum_up=STIndsum_up+STInd_up; end 
            if nout > 2, STIndsum_down=STIndsum_down+STInd_down; end
            if nout > 3, VTIndsum_up=VTIndsum_up+VTInd_up; end 
            if nout > 4, VTIndsum_down=VTIndsum_down+VTInd_down; end 
        end  % end for N=1:Nsim,
        
        TTInd_upp=TTIndsum_upp/Nsim;% max(1,TTIndsum_upp(end));
        TTF_up(kk,:)=TTInd_upp;  
        
        if nout > 0, 
            TTInd_downp=TTIndsum_downp/Nsim; 
            TTF_do(kk,:)=TTInd_downp;        
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
    
    if nout > 0, 
        TTFDO=simpson(uu',wxpdf.*TTF_do);
        TTFDO=TTFDO/TTFDO(end);  
        TTFd.f{vv}=TTFDO(1:end-1);  
    end
    if nout > 1, 
        STFUP=simpson(uu',wxpdf.*STF_up);
        STFUP=STFUP/STFUP(end);
        STFu.f{vv}=STFUP(1:end-1);
    end
    if nout > 2,
        STFDO=simpson(uu',wxpdf.*STF_do);
        STFDO=STFDO/STFDO(end);    
        STFd.f{vv}=STFDO(1:end-1);
    end
    if nout > 3,
        VTFUP=simpson(uu',wxpdf.*VTF_up);
        VTFUP=VTFUP/VTFUP(end);
        VTFu.f{vv}=VTFUP(1:end-1);
    end
    if nout > 4,
        VTFDO=simpson(uu',wxpdf.*VTF_do);
        VTFDO=VTFDO/VTFDO(end); 
        VTFd.f{vv}=VTFDO(1:end-1);
    end
end
if nout > 0, varargout(1)={TTFd}; end
if nout > 1, varargout(2)={ST Fu}; end
if nout > 2, varargout(3)={STFd}; end
if nout > 3, varargout(4)={VTFu}; end
if nout > 4, varargout(5)={VTFd}; end
end

function [TTInd_upp,TTInd_upm,TTInd_downp,TTInd_downm,...
          STInd_up,STInd_down, ...
          VTInd_up,VTInd_down]=g_y(v,yyy,u,SYZ,SZZ,CH,NC,nout)
    noll=[0 0 0 0];
    mu=[0 0 0 1]' + SYZ*(SZZ\[v -u]');
    mu=mu';
    Ident=eye(4);     
    
    U=rndnormnd(noll,Ident,NC);
    Y=U;
    
%    U=mvnrnd(noll,Ident,NC);
%    Y=U';
    
    Y=(CH'*Y)';
    Y=ones(NC,1)*mu+Y;
    Yp=Y(Y(:,4)>0,:); 
    Ym=Y(Y(:,4)<0,:); 

    [NCp,~]=size(Yp); % Exclude retrograd x
    [NCm,~]=size(Ym); % Exclude retrograd x

    TTder=Y(:,1)-Y(:,2).*Y(:,3)./Y(:,4);
    STder=Y(:,2)./Y(:,4);
    VTder=Y(:,3);
    det=abs(Y(:,1).*Y(:,4)-Y(:,2).*Y(:,3));
    Nyyy=length(yyy);
    dety=det*ones(1,Nyyy);

    TTdery_up=(TTder*ones(1,Nyyy))<(ones(NC,1)*yyy);
    TTdery_up=TTdery_up.*(TTder*ones(1,Nyyy)>0);
    dTT=dety.*TTdery_up;
%    TTInd_upp=mean(dTT);
    if NCp>0,
        TTInd_upp=mean(dTT(Y(:,4)>0,:));
    else
        TTInd_upp=zeros(1,Nyyy);
    end
    TTInd_upm=sum(dTT(Y(:,4)<0,:)); 
    
    if nout > 0,
        TTdery_do=(TTder*ones(1,Nyyy))>(-ones(NC,1)*yyy);
        TTdery_do=TTdery_do.*(TTder*ones(1,Nyyy)<0);
        dTT=dety.*TTdery_do;
        if NCp>0,
            TTInd_downp=mean(dTT(Y(:,4)>0,:));
        else
            TTInd_downp=zeros(1,Nyyy);
        end
        TTInd_downm=sum(dTT(Y(:,4)<0,:));
    else
        TTInd_downp=[];
        TTInd_downm=[];
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

