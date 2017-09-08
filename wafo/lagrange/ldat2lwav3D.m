function L=ldat2lwav3D(W,X,Y,opt3D,varargin)
%LDAT2LWAV3D Generates Lagrange 3D wave process from simulated components
%   from W,X,Y fields (as output from spec2ldat3D) 
% 
%CALL:  L =ldat2lwav3D(W,X,Y,opt3D)
%
%   L       = Lagrange structure with fields 
%               L.Z, L.x, L.y, L.t, L.type
%
%   W       = Gaussian vertical process structure w.Z, w.u, w.t
%   X/Y     = Gaussian horizontal variation structures 
%             X.Z, X.u, X.v, X.t
%             Y.Z, Y.u, Y.v, Y.t
%             if  X and Y  are empty then output is Gaussian
%   opt3D   = structure (set by genoptset) with fields 
%       .type   = 'movie' gives time dependent wave fields over times W.t
%               = 'field' give wave field(s) at time(s) given by
%       .t0     't0' (string) or empty: gives t0 = W.t(end)/2 (default)
%               = t0 (numeric): gives field at time  t0
%       .start  = [startx starty]  lower left corner of fields 
%       .end    = [endx endy]  upper right corner of fields 
%                   if empty  end  is generated from  start  coordinates
%       .plotflag   = 'on'  plots one field - the last one
%
%               OR
%       .type   = 'timeseries' gives  n  time series at points
%                    with coordinates
%       .PP     = 2 x n  array [p1,...,pn; q1,...,qn]'
%                   If  n > 1 then output  L.Z  is is 
%                   a cell array  L.Z{1}, ..., L.Z{n} 
%       .rate       = interpolation rate (default = 1) not yet implemented
%       .plotflag   = 'on'  plots one time series - the last one
%
%               OR (NOT YET AVAILABLE)
%       .type   = 'swath' gives encountered Lagrange sea elevation
%                    observed from a moving object with speed(s) --
%                    possible further fields for this option are 
%       .v [m/s]= along straight line(s) from 
%       .start  = default = (W.u(1),W.v(1)/2)  to 
%       .end      default = (W.u(end),W.v(1)/2)
%               OR (NOT YET AVAILABLE)
%       .type   = 'vfield' gives velocity field         
%

%  Tested on Matlab 8.1, 8.6
%  History
%  Created by GL Oct 2014
%  Modified help text and options March 2015
%  Replaced NaN's in field and movie by nearest valid point value, October 2015

if nargin>4,
    opt3D=genoptset(opt3D,varargin{:});
end

fs=14; %FontSize
L=[];
if isempty(X) && isempty(Y) % Only Gaussian model
    model='Gaussian model';
    w=W;
    x.Z=permute(reshape(repmat(W.u',1,length(W.v)*length(W.t)),...
        length(W.u),length(W.v),length(W.t)),[2 1 3]);
    x.u=W.u;
    x.v=W.v;
    x.t=W.t; 
    y.Z=reshape(repmat(W.v,length(W.u)*length(W.t),1),...
        length(W.v),length(W.u),length(W.t));
    y.u=W.u;
    y.v=W.v;
    y.t=W.t;
else % Lagrange modewl
    model='Lagrange model';
    w=W;
    xP=permute(reshape(repmat(X.u',1,length(X.v)*length(X.t)),...
        length(X.u),length(X.v),length(X.t)),[2 1 3]);
    x.Z=X.Z+xP;
    x.u=X.u;
    x.v=X.v;
    x.t=X.t;

    yP=reshape(repmat(Y.v,length(Y.u)*length(Y.t),1),...
        length(Y.v),length(Y.u),length(Y.t));
    y.Z=Y.Z+yP;    
    y.u=Y.u;
    y.v=Y.v;
    y.t=Y.t;
end

if strcmpi(opt3D.type(1:4),'time'), % Gives time series at points  P1,...,Pn
    % First define interpolation rate
    tic
    if ~isfield(opt3D,'rate') || isempty(opt3D.rate)
        rate=1;
    else
        rate=max(1,floor(opt3D.rate));
    end
    % End define interpolation rate

    if ~isfield(opt3D,'PP') || isempty(opt3D.PP), 
        PP = [W.u(end)/2;W.v(end)/2];
    else    
        PP = opt3D.PP;
    end
    [mp,np] = size(PP); 
    if mp==1 && np==2,
        PP = PP';
    elseif mp>=2,
        PP=PP(1:2,:); % Points shall be defined as an  2 x np array
    end
    
    L = struct('t',w.t','Z',[],'P',PP);
    L.Z = cell(1,np);

    [xI,yI]=meshgrid(PP(1,:),PP(2,:));
    nxy=numel(w.Z(:,:,length(w.t)));
 WT=waitbar(0,'Looping over  t');   
    for k=1:length(w.t),
        waitbar(k/length(w.t),WT)
        xx=reshape(squeeze(x.Z(:,:,k)),nxy,1);
        yy=reshape(squeeze(y.Z(:,:,k)),nxy,1);
        ww=reshape(squeeze(w.Z(:,:,k)),nxy,1);
        Vq=griddata(xx,yy,ww,xI,yI);
        for p=1:np,
            L.Z{p}(k)=Vq(p,p);
        end
    end 
 close(WT)
    toc

    if strcmp(opt3D.plotflag,'on'),
        figure(1)
        plot(L.t,L.Z{end})
    end
    if np==1,
        L=cell2struct(L.Z,'Z');
        L.t=w.t;
        L.P=PP;
    end 
end

if strcmpi(opt3D.type(1:5),'field') 
 % Define time and area
    tic
    if ~isfield(opt3D,'t0') || isempty(opt3D.t0)
        t0_index=floor(length(W.t)/2);
        t0=W.t(end)/2;
    elseif isnumeric(opt3D.t0)
        t0=opt3D.t0;
        t0_index=find(W.t>=t0,1);
    end
elseif strcmpi(opt3D.type(1:5),'movie'), 
    t0=W.t;
    t0_index=1:length(t0);
end

if strcmpi(opt3D.type(1:5),'field') || strcmpi(opt3D.type(1:5),'movie'), 
    tic
    if ~isfield(opt3D,'start') || isempty(opt3D.start)
        start=[W.u(5) W.v(5)];
        startindex=[5 5];
    else
        start=opt3D.start;
        startindex=[find(W.u>=opt3D.start(1),1) find(W.v>=opt3D.start(2),1)];
    end
    if ~isfield(opt3D,'end') || isempty(opt3D.end)
        slut=[W.u(end)-start(1) W.v(end)-start(2)];
        slutindex=[length(W.u)-startindex(1)+1  length(W.v)-startindex(2)+1];
    else
        slut=opt3D.end;
        slutindex=[find(W.u>=opt3D.end(1),1) find(W.v>=opt3D.end(2),1)];
    end
 % End Define time and area
   %startindex
   %slutindex
    L.x=W.v(startindex(2):slutindex(2)); %W.v;
    L.y=W.u(startindex(1):slutindex(1)); %W.u;
    L.t=t0;
    [xI,yI]=meshgrid(W.u(startindex(1):slutindex(1)),W.v(startindex(2):slutindex(2)));
    nxy=numel(w.Z(:,:,1));

    L.Z=zeros(length(xI(:,1)),length(yI(1,:)),length(t0));
    WT=waitbar(0,'Looping over t');

    for t=1:length(t0),
    waitbar(t/length(t0),WT)
        xx=reshape(squeeze(x.Z(:,:,t0_index(t))),nxy,1);
        yy=reshape(squeeze(y.Z(:,:,t0_index(t))),nxy,1);
        ww=reshape(squeeze(w.Z(:,:,t0_index(t))),nxy,1);
        Vq=griddata(xx,yy,ww,xI,yI);
        
        if sum(find(isnan(Vq)))>0,
            disp('Simulated field has moved outside study area')
            disp('Interpolation inaccurate or sometimes impossible')
            [I1,I2] = ind2sub(size(Vq),find(isnan(Vq)));  
            [J1,J2] = ind2sub(size(Vq),find(~isnan(Vq(:,1:2))));
            IDX = dsearchn([J1 J2], [I1 I2]); % Find nearest valid point
            for miss=1:length(I1), % Replace NaN's
                Vq(I1(miss),I2(miss))=Vq(J1(IDX(miss)),J2(IDX(miss)));
            end
        end
        L.Z(:,:,t) = Vq;
    end

    close(WT)
    if strcmp(opt3D.plotflag,'on'), % Plotting last field
        figNo = gcf;
        figure(figNo)  
        colormap('default')
        surf(xI,yI,5*L.Z(:,:,end))
        colorbar
        axis equal
        shading interp
        set(gca,'FontSize',fs)
        zlabel('5 \times L.Z','FontSize',fs-1)
    end
    
    toc
end
if strcmpi(opt3D.type(1:3),'swa'), % Not yet implemented
    v=opt3D.v;
    if isempty(opt3D.start),
    end
end
return
