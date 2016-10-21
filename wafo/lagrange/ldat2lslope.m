function Slopes=ldat2lslope(w,x,typ,levels)
%LDAT2LSLOPE Extracts slopes at level crossings in Lagrange model
%
%Call: Slopes  = ldat2lslope(w,x,type,levels)
%               
%   Slopes  = struct array with observed slopes at the up- and 
%             downcrossings of specified levels
%
%   w,x     = vertical and horizontal component in Lagrange model
%   type    = 'space' gives slopes in space waves
%                 'time' gives time slopes
%   levels  = vector of absolute levels relative mwl=0 
%                 (no default)
%
% Example:
%   S=jonswap; mom=spec2mom(S);
%   opt=simoptset;
%   opt=simoptset(opt,'dt',0.25,'du',0.25)
%   [w,x]=spec2ldat(S,opt)
%   levels=[0 1 2]*sqrt(mom(1));
%   Slopes=ldat2lslope(w,x,'time',levels)

% Tested on Matlab 8.1, 8.6
%   History:
%   Modified for time 2015-03-02 to better identify x.Z/w.Z crossings
%   Updated 2015-02-02 to Matlab standards by GL
%   New time version 2014-02-26 by GL
%   Included time waves 2008-11-15, GL
%   Created 2007-012-16 for space waves by GL

if nargin<4,
    disp('Absolute levels required')
    return
end

nlevel=length(levels);
Slopes = struct('up',[],'down',[],'levels',levels,'type',typ);
for ind=1:nlevel,
    Slopes.up{ind}=[];
    Slopes.down{ind}=[];
end

X=x.u*ones(size(x.Z(1,:)))+x.Z;
W=w.Z;
if strcmp(typ,'time')
    u0=floor((w.u(1)+w.u(end))/2);
    dt=w.t(2)-w.t(1);
    du=w.u(2)-w.u(1);
    C=contourc(x.t,x.u,X,[u0 u0]);
    Cx=C(1,2:C(2,1)+1);
    Cy=C(2,2:C(2,1)+1);
    Cxi=[fix(Cx/dt); Cx/dt-fix(Cx/dt)];
    Cyi=[fix(Cy/du); Cy/du-fix(Cy/du)];
    Cti=min([Cxi(1,1:end-1);Cxi(1,2:end)]); 
    Cui=min([Cyi(1,1:end-1);Cyi(1,2:end)]); 
    Cti=Cti+1; % Lower left corner t-index of x.Z-crossing rectangle 
    Cui=Cui+1; % Lower left corner u-index of x.Z-crossing rectangle
    Nruta=length(Cti);
    Xval=zeros(4,Nruta);
    Wval=zeros(4,Nruta);
    
    for kol=1:Nruta, % Get x- and w-values at rectangle corners
        Xval(:,kol)=...
            [X(Cui(kol),Cti(kol)) X(Cui(kol)+1,Cti(kol))...
             X(Cui(kol)+1,Cti(kol)+1) X(Cui(kol),Cti(kol)+1)]';
        Wval(:,kol)=...
            [W(Cui(kol),Cti(kol)) W(Cui(kol)+1,Cti(kol))...
             W(Cui(kol)+1,Cti(kol)+1) W(Cui(kol),Cti(kol)+1)]';
    end
    
    for v=1:nlevel,
        lv=levels(v);
        Possible=(max(Wval(:,:))>lv).*(min(Wval(:,:))<lv); 
            % Selects rectangles with w.Z-crossings 
        Impossible=... % Find saddle points
            (sign(Wval(1,:)-lv)==sign(Wval(3,:)-lv)).*...
            (sign(Wval(2,:)-lv)==sign(Wval(4,:)-lv));
        Possible=Possible.*(1-Impossible); % Remove saddle points 
            % Next: Select rectangle where x- and w-curves cross
        Ctip=Cti(Possible==1); % Lower left corner t-index of possible rectangle 
        Cuip=Cui(Possible==1); % Lower left corner u-index of possible rectangle 
        Npossible=length(Ctip);
        RektI=find(Possible);

        for r=1:Npossible,
            Crw=contourc([x.t(Ctip(r)) x.t(Ctip(r)+1)],[x.u(Cuip(r)) x.u(Cuip(r)+1)]',...
                [W(Cuip(r),Ctip(r)) W(Cuip(r),Ctip(r)+1);...
                 W(Cuip(r)+1,Ctip(r)) W(Cuip(r)+1,Ctip(r)+1)],[lv lv]); 
            Crx=[Cx(RektI(r)) Cx(RektI(r)+1); Cy(RektI(r)) Cy(RektI(r)+1)];
            Tp=(Crw(2,2)-Crx(2,1)-(Crw(1,2)-Crx(1,1))*...
                (Crx(2,1)-Crx(2,2))/(Crx(1,1)-Crx(1,2)));
            Tq=(Crw(2,3)-Crx(2,1)-(Crw(1,3)-Crx(1,1))*...
                (Crx(2,1)-Crx(2,2))/(Crx(1,1)-Crx(1,2)));
            if Tp>0 && Tq<0 || Tp<0 && Tq>0
%               do nothing 
            else
                Possible(RektI(r))=0;
            end
        end
        Xcross=Xval(:,Possible==1);
        Wcross=Wval(:,Possible==1);
        Wt=(Wcross(3,:)-Wcross(2,:)+Wcross(4,:)-Wcross(1,:))/2/dt;
        Xt=(Xcross(3,:)-Xcross(2,:)+Xcross(4,:)-Xcross(1,:))/2/dt;
        Wu=(Wcross(3,:)-Wcross(4,:)+Wcross(2,:)-Wcross(1,:))/2/du;
        Xu=(Xcross(3,:)-Xcross(4,:)+Xcross(2,:)-Xcross(1,:))/2/du;
        Slop=Wt-Wu.*Xt./Xu;
        Slop=Slop(Xu>0); % Exclude retrograd slopes
        Slopes.up{v}=Slop(Slop>0);
        Slopes.down{v}=Slop(Slop<0);          
    end
    Slopes.note='Slopes in time waves';   
    Slopes.levels=levels;
end

if strcmp(typ,'space'),
    t0=floor((w.t(1)+w.t(end))/2);
    W=w.Z(:,t0);
    X=x.Z(:,t0);
    dd=w.u(2)-w.u(1);

    for j=1:nlevel,
        [ind_up,~]=dat2crossind(W,levels(j),'u',true);
        [ind_down,~]=dat2crossind(W,levels(j),'d',true);
        Wu=gradient(W,dd);
        Xu=1+gradient(X,dd);
        Wd=gradient(W,dd);
        Xd=1+gradient(X,dd);
        Sup=Wu./Xu;
        Sdo=Wd./Xd;
        Slopes.up{j}=max(0,Sup(ind_up));
        Slopes.down{j}=min(0,Sdo(ind_down));
    end
    Slopes.note='Slopes in space waves';
    Slopes.levels=levels;
end
return

