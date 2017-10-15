% Script with commands for WafoL tutorial, Chapter 3

rng('default')
pstate = pause('off');
TSTART = cputime;
disp(' Chapter 3: 2nd order non-linear Lagrange waves')

disp(' Section 3.2.1')
    % Figure: Directional spectrum
    S = jonswap;
    D = spreading(linspace(-pi,pi,91),'cos2s',0,10,S.w,0);
    S.h = 40; 
    S.D = D;
    Sdir = mkdspec(S,D);
    figure(1), clf
    plotspec(Sdir)
    option = simoptset('Nt',1001,'dt',0.2,'Nu',300,'du',1,'Nv',100,'dv',1);
    pause
    
disp(' Section 3.2.2')
    % Figure: qq-plot
    [W1,X1,Y1] = spec2ldat3D(Sdir,option);
    [W2,X2,Y2] = spec2ldat3DM(S,1,option);

    % Run if PCT is available
    % pool=parpool(4);
    % [W3,X3,Y3] = spec2ldat3DP(S,1,option);

    figure(2), clf
    subplot(221)
    plotqq(W1.Z(:),W2.Z(:)); grid on
    xlabel('W1.Z quantiles')
    ylabel('W2.Z quantiles')

    subplot(222)
    plotnorm(W1.Z(:))
    pause
    clear W1 W2 X1 X2 Y1 Y2
    
disp(' Section 3.2.3')
    % Check changes in region
    S = jonswap;
    opt=simoptset('Nt',501,'dt',0.5,'Nu',501,'du',1);
    [W,X] = spec2ldat(S,opt);
    [Wm,Xm] = spec2ldat3DM(S,1,opt);
    pause

disp(' Section 3.3.1')
    % Figure: Linear and non-lines with spec2nlsdat
    S = jonswap; S.h = 20;
    np = 250; dt = 0.2;
    [xs2, xs1] = spec2nlsdat(S,np,dt);
    figure(1), clf
    waveplot(xs1,'b',xs2,'r',1,1)
    pause

disp(' Section 3.3.2')
    % Time waves
    S = jonswap; S.h = 20;
    opt = simoptset('Nt',250,'dt',0.2,'Nu',1001,'du',1,'Nv',1);
    [W,X,~,W2,X2,~] = spec2ldat3DM(S,2,opt);
    X2 % The 2nd order x-components has a field 'drift'
    pause
    
    % Figure: Time waves with and without Stokes drift
    Wtot2 = W; Wtot2.Z = W.Z + W2.Z;
    Wtot2d = Wtot2;  % No Stokes drift in the vertical direction ! 
    Xtot2 = X; Xtot2.Z = X.Z + X2.Z;
    Xtot2d = X; Xtot2d.Z = X.Z + X2.Z + X2.drift;

    [L,L0] = ldat2lwav(W,X,'time',250,1,0);
    [L2,L20] = ldat2lwav(Wtot2,Xtot2,'time',250,1,0);
    [L2d,L20d] = ldat2lwav(Wtot2d,Xtot2d,'time',250,1,0);

    figure(1), clf
    plot(L0.t,L0.Z,'b','LineWidth',2); hold on, grid on;
    plot(L20.t,L20.Z,'r-.','LineWidth',2)
    plot(L20d.t,L20d.Z,'r','LineWidth',2)
    set(gca,'FontSize',14)
    title('Surface elevation from mean water level (MWL)')
    xlabel('Time (sec)'); hold off
    pause

    % Space waves
    S = jonswap; S.h = 20;
    opt = simoptset('Nt',1000,'dt',0.2,'Nu',1001,'du',1,'Nv',1);
    [W,X,~,W2,X2,~] = spec2ldat3DM(S,2,opt);
    
    % Figure: Space wave with Stokes drift
    Wtot2 = W; Wtot2.Z = W.Z + W2.Z;
    Wtot2d = Wtot2;  % No Stokes drift in the vertical direction ! 
    Xtot2 = X; Xtot2.Z = X.Z + X2.Z;
    Xtot2d = X; Xtot2d.Z = X.Z + X2.Z + X2.drift;

    [L50,L050] = ldat2lwav(W,X,'space',50,1,0);
    %[L250,L2050] = ldat2lwav(Wtot2,Xtot2,'space',50,1,0);
    [L2d50,L20d50] = ldat2lwav(Wtot2d,Xtot2d,'space',50,1,0);
    [L100,L0100] = ldat2lwav(W,X,'space',100,1,0);
    %[L2100,L20100] = ldat2lwav(Wtot2,Xtot2,'space',100,1,0);
    [L2d100,L20d100] = ldat2lwav(Wtot2d,Xtot2d,'space',100,1,0);
    [L150,L0150] = ldat2lwav(W,X,'space',150,1,0);
    %[L2150,L20150] = ldat2lwav(Wtot2,Xtot2,'space',150,1,0);
    [L2d150,L20d150] = ldat2lwav(Wtot2d,Xtot2d,'space',150,1,0);

    figure(1), clf
    subplot(311)
    plot(L050.u,L050.Z,'b','LineWidth',2); hold on, grid on;
    plot(L20d50.u,L20d50.Z,'r','LineWidth',2)
    set(gca,'FontSize',12)
    title('Surface elevation from mean water level (MWL)')
    xlabel('Location (m) at time = 50 (sec)')
    subplot(312)
    plot(L0100.u,L0100.Z,'b','LineWidth',2); hold on, grid on;
    plot(L20d100.u,L20d100.Z,'r','LineWidth',2)
    set(gca,'FontSize',12)
    xlabel('Location (m) at time = 100 (sec)')
    subplot(313)
    plot(L0150.u,L0150.Z,'b','LineWidth',2); hold on, grid on;
    plot(L20d150.u,L20d150.Z,'r','LineWidth',2)
    set(gca,'FontSize',12)
    xlabel('Location (m) at time = 150 (sec)')
    pause
       
disp(' Section 3.3.3, Front-back asymmetry')
    % Figure: Front-back asymmetric 2nd order waves
    S = jonswap; S.h = 20;   
    [W,X,~,W2,X2,~] = spec2ldat3DM(S,2,opt,'Nt',250,'lalpha',1.5);
    Wtot2 = W; Wtot2.Z = W.Z + W2.Z;
    Wtot2d = Wtot2;  
    Xtot2 = X; Xtot2.Z = X.Z + X2.Z;
    Xtot2d = X; Xtot2d.Z = X.Z + X2.Z + X2.drift;

    [L,L0] = ldat2lwav(W,X,'time',250,1,0);
    [L2,L20] = ldat2lwav(Wtot2,Xtot2,'time',250,1,0);
    [L2d,L20d] = ldat2lwav(Wtot2d,Xtot2d,'time',250,1,0);

    figure(1), clf
    subplot(311)
    plot(L0.t,L0.Z,'b','LineWidth',2); hold on, grid on;
    plot(L20.t,L20.Z,'r-.','LineWidth',2)
    plot(L20d.t,L20d.Z,'r','LineWidth',2)
    set(gca,'FontSize',14)
    title('Surface elevation from mean water level (MWL)')
    xlabel('Time (sec)')
    pause

disp(' Section 3.4')
    S = jonswap; S.h=20;
    D = spreading(linspace(-pi,pi,51),'cos2s',0,15,S.w,0);
    S.D = D;
    opt=simoptset('Nt',250,'dt',0.2,'Nu',250,'du',1,'Nv',50,'dv',1);
    [W,X,Y,W2,X2,Y2] = spec2ldat3DM(S,2,opt);
    
    % Figure: Check the Stokes drift
    figure(10), clf
    surf(X2.u(1:2:end),X2.v(1:2:end),X2.drift(1:2:end,1:2:end,end)) 
    
    % Extend the reference region
    opt=simoptset('Nt',250,'dt',0.2,'Nu',360,'du',1,'Nv',130,'dv',1);
    [W,X,Y,W2,X2,Y2] = spec2ldat3DM(S,2,opt);
    Wtot2 = W; Wtot2.Z = W.Z + W2.Z;
    Xtot2d = X; Xtot2d.Z = X.Z + X2.Z + X2.drift;
    Ytot2d = Y; Ytot2d.Z = Y.Z + Y2.Z + Y2.drift;
    % and set the region of interest
    opt3D = genoptset('start',[10 15],'end',[260 115]);
    
    % Figures: Fields and movies
    figure(1), clf
    L=ldat2lwav3D(Wtot2,Xtot2d,Ytot2d,opt3D)
    drawnow
    figure(2), clf
    Mv=seamovie(L,1)
    pause

    % Front-back asymmetry
    optasym = simoptset(opt,'lalpha',1.5);
    [W,X,Y,W2,X2,Y2] = spec2ldat3DM(S,2,optasym);
    Wtot2 = W; Wtot2.Z = W.Z + W2.Z;
    Xtot2d = X; Xtot2d.Z = X.Z + X2.Z + X2.drift;
    Ytot2d = Y; Ytot2d.Z = Y.Z + Y2.Z + Y2.drift;
    figure(3), clf
    L2=ldat2lwav3D(Wtot2,Xtot2d,Ytot2d,opt3D)
    drawnow
    figure(4), clf
    Mv2=seamovie(L2,1)

    disp('End of Chapter 3')
       
    TSTOP=cputime
    disp(['Total time = ' num2str(TSTOP-TSTART) ' sec']);
pause(pstate)   
