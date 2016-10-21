function [L,L0]=ldat2lwav(w_in,x_in,type,tu0,dense,plott)
%LDAT2LWAV Finds time/space Lagrange process from simulated components
%   This version returns true time or 2D space profile and  
%   the smoothed initial part until the first loop/break
%                 
%CALL: [L,L0]=ldat2lwav(w,x,type,tu0,dense,plotting)
%
%   L        = Lagrange process structure with fields 
%                  L.type and L.t/L.u and L.Z
%                  L can contain loops and breaking waves
%   L0       = a pruned and smoothed version of L without loops
%
%   w        = Gaussian vertical process structure w.Z, w.u, w.t
%   x        = Gaussian horizontal variation process structure 
%              x.Z, x.u, x.t
%   type     = 'time' or 'space' 
%   tu0      = time t0 for space wave, 
%            = space coordinate u0 for time wave
%   dense    = interpolation rate for smoothing
%   plotting = 0/false, no plotting (default), 
%            = 1/true, plot L (and L0 if not empty)

% Tested on Matlab 8.1, 8.6
% Revised December 2015 to estimate mean period or wave length from data
% Completely new version 2015-02-17 by GL
% Created September 2007 by GL with numerous modifications

L=[];
L0=[];

if nargin<6 || isempty(plott),
    plotting=false;
else 
    plotting=logical(plott);
end

if nargin<5 || isempty(dense),
    dense=10;
end

fs=15;
x.Z=x_in.u*ones(size(x_in.Z(1,:)))+x_in.Z;
x.u=x_in.u;
x.t=x_in.t;
w=w_in;

if strcmp(type,'time'), % Compute exact time waves
    if nargin<4 || isempty(tu0)
        tu0=w.u(end)/2;
    end
    dt=w.t(2)-w.t(1);
    totT=w.t(end)-w.t(1);
    dtt=dt/dense;
    if isrow(w.t); tcol=w.t'; else tcol=w.t; end 
    meanperiod=[];
    for ku=1:length(w.u);
        TC = dat2tc([tcol w.Z(ku,:)'],0);
        meanp=totT/max(1,sum(TC(:,2)>0));
        meanperiod=[meanperiod meanp];
    end
    meanperiod=mean(meanperiod);
    xz_=x.Z(1:end-1,:);
    xz__=x.Z(2:end,:);
    s= (xz_<tu0).*(xz__>=tu0) + (xz_>=tu0).*(xz__<tu0);
    [I,J]=find(s);
    ni=length(I);
    for k=1:ni;
        L.t(k)=x.t(J(k)); 
        L.Z(k)=w.Z(I(k),J(k));
    end
    if plotting
        figure(10)
        subplot(211)
        plot(L.t,L.Z,'.','MarkerSize',5)
        set(gca,'FontSize',fs)
        V=axis;
        grid on
    end
    WW_.t=w.t(1):dtt:w.t(end);
    
% Check if there are double crossings
% and if so, delete the part from the first one
    fd=find(diff(J)==0,1); % First time column with double crossings
    if isempty(fd), % No double crossings
        WW.t=L.t;
        WW.Z=L.Z;
        WW_.t=WW_.t(WW_.t>=min(L.t));
        WW_.t=WW_.t(WW_.t<=max(L.t));
%        WW_.t=(0:dense*length(WW.t))/dense/length(WW.t)*max(WW.t)+min(WW.t);
        WW_.Z=cssmooth(WW.t,WW.Z,1,WW_.t);    
        L0=WW_;
    else
        fd_=fd-round(w.meanperiod/(w.t(2)-w.t(1))/5);
        if logical((fd_>2)) && logical((L.t(fd_) >= 4*meanperiod)),
            WW.t=L.t(1:fd_); 
            WW.Z=L.Z(1:fd_);
            WW_.t=WW_.t(WW_.t>=L.t(1));
            WW_.t=WW_.t(WW_.t<=L.t(fd_));
%            WW_.t=(0:dense*length(WW.t))/dense/length(WW.t)*max(WW.t)+min(WW.t);
            WW_.Z=cssmooth(WW.t,WW.Z,1,WW_.t);    
            L0=WW_;
        else
            L0=[];
        end
    end
    if plotting,
        subplot(212)
        if ~isempty(L0)
            plot(L0.t,L0.Z,'LineWidth',1.5)
            set(gca,'FontSize',fs) 
            axis(V)
            grid on
        end
    end

elseif strcmp(type,'space'), % Space waves
    if nargin<4|isempty(tu0)
        tu0=w.t(end)/2;
    end
    du=w.u(2)-w.u(1);
    totU=w.u(end)-w.u(1);
    duu=du/dense;
    meanwavelength=[];
    for kt=1:length(w.t);
        TC = dat2tc([w.u w.Z(:,kt)],0);
        meanl=totU/max(1,sum(TC(:,2)>0));
        meanwavelength=[meanwavelength meanl];
    end
    meanwavelength=mean(meanwavelength);
%   Shift location to make space wave at time  tu0
    J=find((w.t>tu0),1)-1; 
    if isempty(J), 
        J=length(w.t);
    else
        J=max(1,J);
    end
    
    WW.Z=w.Z(:,J);
    x0.Z=x.Z(:,J);
    WW.u=x0.Z;
    L=WW;
    if plotting,
        figure(10)
        subplot(211)
        plot(L.u,L.Z,'LineWidth',1.5); 
        set(gca,'FontSize',fs)
        V=axis;
        grid on
    end
    WW_.u=w.u(1):duu:w.u(end);
    
% Check if there are double crossings
% and if so, delete the part from the first one
    fd=find(diff(L.u)<0,1); % First space row with backward movement
    if isempty(fd),
        WW.u=L.u;
        WW.Z=L.Z;
        WW_.u=WW_.u(WW_.u>=min(L.u));
        WW_.u=WW_.u(WW_.u<=max(L.u));
%        WW_.u=(0:dense*length(WW.u))/dense/length(WW.u)*max(WW.u)+min(WW.u);
        WW_.Z=cssmooth(WW.u,WW.Z,1,WW_.u); 
        WW_.u=WW_.u';
        WW_.Z=WW_.Z';
        L0=WW_;
    else
        fd_=fd-round(w.meanwavelength/(w.u(2)-w.u(1))*2/3);
        if logical((fd_>2)) && logical((L.u(fd_) >= 4*meanwavelength)),
            WW.u=L.u(1:fd_); 
            WW.Z=L.Z(1:fd_);
            WW_.u=WW_.u(WW_.u>=min(WW.u));
            WW_.u=WW_.u(WW_.u<=max(WW.u));
%            WW_.u=(0:dense*length(WW.u))/dense/length(WW.u)*max(WW.u)+min(WW.u);
            WW_.Z=cssmooth(WW.u,WW.Z,1,WW_.u);  
            WW_.u=WW_.u';
            WW_.Z=WW_.Z';
            L0=WW_;
        else
            L0=[];
        end
    end
    if plotting,
        subplot(212)
        if ~isempty(L0)
            plot(L0.u,L0.Z,'LineWidth',1.5); 
            set(gca,'FontSize',fs)
            axis(V);
            grid on
        end
    end
else 
    disp('Unknown wave type')
end
return

