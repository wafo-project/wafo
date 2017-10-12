% Script with commands for WafoL tutorial, Chapter 2

rng('default')
pstate = pause('off');
TSTART = cputime;
disp(' Chapter 2: 3D Lagrange waves')

disp(' Section 2.2.1  Spectrum and simulation options')
   % Directional spectrum
   % Figure: Directional spectra in polar and Cartesian form
   S = jonswap
   D = spreading(linspace(-pi,pi,51),'cos2s',0,15,S.w,0)
   figure(1); clf 
   Snew = mkdspec(S,D,1) 

   Sk2d = spec2spec(Snew,'k2d')
   figure(2); clf
   plotspec(Sk2d); clear Sk2d
   axis([-0.05 0.05 -0.05 0.05])
   pause
   
   % Simulation options
   optTime = simoptset('Nt',1001,'dt',0.2,...
                 'Nu',128,'du',5,'Nv',64,'dv',5);
   optSpace = simoptset('Nt',1,'dt',0.5,...
                 'Nu',1024,'du',1,'Nv',512,'dv',1);

   % Generation of the elementary processes
   S = jonswap(1.5); S.h = 20;
   D5 = spreading(101,'cos',0,5,S.w,0);
   SD5 = mkdspec(S,D5);
   D15 = spreading(101,'cos',0,15,S.w,0); 
   SD15 = mkdspec(S,D15);

%   [W,X,Y] = spec2ldat3D(SD5,optTime);
   [W,X,Y] = spec2ldat3D(SD5,optSpace);
   pause
   
disp(' Section 2.2.2  Generating the Lagrange waves from the 3D fields')

   % Generating a single field
   optSpace = simoptset('Nt',20,'dt',1,...
                          'Nu',256,'du',1,'Nv',128,'dv',1);
   opt3D = genoptset('type','field','t0',10)

   % Figure: One front-back asymmetric field
   [W,X,Y] = spec2ldat3D(SD15,optSpace,'lalpha',1)
   Sx = mean(mean(std(X.Z))) % result = 4.6
   Sy = mean(mean(std(Y.Z))) % result = 1.9
   opt3D = genoptset(opt3D,'start',[20 10])
   Lfield = ldat2lwav3D(W,X,Y,opt3D)
   pause
   
   % Generating a movie
   optTime = simoptset('Nt',101,'dt',0.2,...
                          'Nu',128,'du',10,'Nv',64,'dv',10);
   opt3D = genoptset('type','movie')

   [W,X,Y] = spec2ldat3D(SD15,optTime,'lalpha',1.5)
   opt3D = genoptset(opt3D,'start',[20 10])
   Lmovie = ldat2lwav3D(W,X,Y,opt3D)
   pause
   
   % Figure: Last frame of seamovie - two displays
   Mv2 = seamovie(Lmovie,2)
   pause
   Mv3 = seamovie(Lmovie,3)
   pause
	
   % Generating time series
   opt3D=genoptset('type','timeseries','PP',[400 425 450; 300 300 300])

   % Figure: Three point time series at distance
   Lseries = ldat2lwav3D(W,X,Y,opt3D,'plotflag','off')
   subplot(311)
   plot(Lseries.t,Lseries.Z{1},'r')
   grid on; hold on
   plot(Lseries.t,Lseries.Z{2},'g')
   plot(Lseries.t,Lseries.Z{3},'b')
   pause
   
   % Figure: Cross-correlation function between time series
   Rx = xcov(Lseries.Z{1},Lseries.Z{3},200,'biased');
   subplot(211)
   plot((-200:200)/5,Rx)
   
   disp('End of Chapter 2')

   TSTOP=cputime
   disp(['Total time = ' num2str(TSTOP-TSTART) ' sec']);
pause(pstate)
