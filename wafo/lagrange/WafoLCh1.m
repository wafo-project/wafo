% Script with commands for WafoL tutorial, Chapter 1

rng('default')
pstate = pause('off');
TSTART = cputime;
disp(' Chapter 1: 2D Lagrange waves')

disp(' Section 1.1.2')
% Figure: Lagrange time wave, example of asymmetry
   S = jonswap; S.h = 20;
   opt = simoptset('dt',1);
   [w,x] = spec2ldat(S,opt,'lalpha',1);
   [L,L0] = ldat2lwav(w,x,'time',[],10);
   subplot(211)
   plot(L.t,L.Z); hold on 
   plot(L0.t,L0.Z,'r')
   axis([0 40 -6 6]); hold off
   pause
   
disp(' Section 1.2.2')
   % Figure: Lagrange space wave, direct plot
   S = jonswap(1.5); S.h = 8;
   opt = simoptset('Nt',256,'dt',0.125,'Nu',...
      256*8,'du',0.25,'iseed',123791);
   [w,x] = spec2ldat(S,opt,'iseed','shuffle')
   subplot(211)
   plot(x.u+x.Z(:,128),w.Z(:,128))
   axis([0 500 -10 10])
   pause
 
disp(' Section 1.2.3')
   % Figure: Lagrange space wave from [L,L0]
   [L,L0] = ldat2lwav(w,x,'space',16)
   subplot(212)
   plot(L.u,L.Z,L0.u,L0.Z,'r')
   axis([0 500 -10 10])
   pause
   
   MGauss = mean(w.Z(:))
   MLagrange = mean(L0.Z)
   SGauss = std(w.Z(:))
   SLagrange = std(L0.Z)
   pause
   
disp(' Section 1.2.4')
   % Figure: Crest-trough asymmetry in space, depth dependence
   opt = simoptset('Nt',128,'Nu',2048,'du',0.25);
   S = jonswap(1.5);
   S3 = S; S3.h = 3;
   S8 = S; S8.h = 8;
   S32 = S; S32.h = 32;
   [w,x] = spec2ldat(S,opt);
   [w3,x3] = spec2ldat(S3,opt);
   [w8,x8] = spec2ldat(S8,opt);
   [w32,x32] = spec2ldat(S32,opt);
   [L,L0] = ldat2lwav(w,x,'space');
   [L3,L03] = ldat2lwav(w3,x3,'space'); 
   [L8,L08] = ldat2lwav(w8,x8,'space');
   [L32,L032] = ldat2lwav(w32,x32,'space');
   figure(1)
   clf
   subplot(411)
   plot(L0.u,L0.Z); axis([0 500 -5 5]);
   title('Depth = \infty')
   subplot(412)
   plot(L032.u,L032.Z); axis([0 500 -5 5]);
   title('Depth = 32 m')
   subplot(413)
   plot(L08.u,L08.Z); axis([0 500 -5 5]); 
   title('Depth = 8 m')
   subplot(414) 
   % 3m water depth may cause loops and ldat2lwav give empty L3,L03
   % We plot space wave directly from w3,x3
   plot(x3.u+x3.Z(:,64),w3.Z(:,64))
   axis([0 500 -5 5])
   title('Depth = 3 m')
   clear w x w3 x3 w8 x8 w32 x32
   pause
   
   % Figure: Truncation of time and space waves  
   S = jonswap(1.5); S.h=20;
   opt = simoptset('dt',0.125,'lalpha',2);
   [w,x] = spec2ldat(S,opt);
   [L,L0] = ldat2lwav(w,x,'time',[],10,1)
   pause
   S = jonswap(1.5); S.h=4;
   opt = simoptset('dt',0.125,'lalpha',0,'ffttype','fftspace');
   [w,x] = spec2ldat(S,opt);
   [L,L0] = ldat2lwav(w,x,'space',[],10,1)
   clear w x 
   pause
   
disp(' Section 1.2.5')
   % Figure: Front-back asymmetric time waves, different alpha
   opt = simoptset('Nt',2048,'dt',0.125,'Nu',512,'du',0.25);
   S = jonswap(1.5); S.h=20;
   [w0,x0] = spec2ldat(S,opt);
   [w1,x1] = spec2ldat(S,opt,'lalpha',0.75);
   [w2,x2] = spec2ldat(S,opt,'lalpha',1.5);
   [L,L0] = ldat2lwav(w0,x0,'time');
   [L1,L01] = ldat2lwav(w1,x1,'time');
   [L2,L02] = ldat2lwav(w2,x2,'time');
   figure(1)
   clf
   subplot(311)
   plot(L0.t,L0.Z)
   title('\alpha = 0')
   axis([0 50 -10 10])
   subplot(312)
   plot(L01.t,L01.Z)
   title('\alpha = 0.75')
   axis([0 50 -10 10])
   subplot(313)
   if ~isempty(L02),
       plot(L02.t,L02.Z)
       title('\alpha = 1.5')
       axis([0 50 -10 10])
   end
   pause
   
disp(' Section 1.3.1')
   % Figure: Empirical time slope CDF at crosings - different alpha 
   opt = simoptset('Nt',2048*16,'dt',0.125,'Nu',256,'du',0.25);
   S = jonswap(1.5,[6 10]); S.h=20;
   [w0,x0] = spec2ldat(S,opt);
   [w1,x1] = spec2ldat(S,opt,'lalpha',0.5);
   [w2,x2] = spec2ldat(S,opt,'lalpha',1);
   [L0,L00] = ldat2lwav(w0,x0,'time'); clear w0 x0
   [L1,L01] = ldat2lwav(w1,x1,'time'); clear w1 x1
   [L2,L02] = ldat2lwav(w2,x2,'time'); clear w2 x2
   mom = spec2mom(S);
   levels=[0 1 2]*sqrt(mom(1)); % wav2slope requires absolute levels
   Slope0 = wav2slope(L00,levels); clear L00
   Slope1 = wav2slope(L01,levels); clear L01
   Slope2 = wav2slope(L02,levels)
   clear L02
   pause
   
   % Plotting
   
   figure(10)
   clf
   plotedf(Slope0.up{1}); hold on
   %get(H101,'Children')
   plotedf(Slope1.up{1});
   plotedf(Slope2.up{1}); 
   plotedf(-Slope0.down{1},'r-.')
   plotedf(-Slope1.down{1},'r-.');
   plotedf(-Slope2.down{1},'r-.'); 
   axis([0 8 0 1])
   pause
   
   % Figure: Empirical space slope CDF at crossings - different alpha
%   S=jonswap(1.5); S.h=20;
   S=jonswap(1.5,[6 10]); S.h=20;
   opt=simoptset('Nt',64,'dt',0.25,...
        'Nu',2048*16,'du',0.25,'ffttype','fftspace');
   Nsim=10;
   Slopes = spec2slopedat(S,Nsim,'space',[],opt)
   Slopes1 = spec2slopedat(S,Nsim,'space',[],opt,'lalpha',1)
   figure(1);   clf
   subplot(221); box;  hold on
   for f=1:length(Slopes.levels),
      if ~isempty(Slopes.up{f}) 
         plotedf(Slopes.up{f}); grid on
      end
      if ~isempty(Slopes.down{f}) 
         plotedf(-Slopes.down{f},'-.'); grid on
      end
   end
   axis([0 0.8 0 1])
   title('\alpha = 0')
   subplot(222); box;  hold on
   for f=1:length(Slopes1.levels),
      if ~isempty(Slopes1.up{f}) 
         plotedf(Slopes1.up{f}); grid on
      end
      if ~isempty(Slopes1.down{f}) 
         plotedf(-Slopes1.down{f},'-.'); grid on
      end
   end
   axis([0 0.8 0 1])
   title('\alpha = 1')
   subplot(223); box;  hold on
   plot(Slopes.meanwavex,Slopes.meanwaveup)
   ax=axis;
   axis([Slopes.meanwavex(1) Slopes.meanwavex(end) ax(3) ax(4)]);
   plot(Slopes.meanwavex,Slopes.meanwavedown,'-.'); grid on 
   subplot(224); box;  hold on
   plot(Slopes1.meanwavex,Slopes1.meanwaveup)
   axis([Slopes.meanwavex(1) Slopes.meanwavex(end) ax(3) ax(4)]);
   plot(Slopes1.meanwavex,Slopes1.meanwavedown,'-.'); grid on
   pause
   
   % Figure: Same as previous, but for time waves
   S=jonswap(1.5,[6 10]); S.h=20;
   opt=simoptset('Nt',2048*8,'dt',0.125,...
        'Nu',321,'du',0.125,'ffttype','ffttime');
   Nsim=10;
   Slopes=spec2slopedat(S,Nsim,'time',[],opt)
   Slopes1=spec2slopedat(S,Nsim,'time',[],opt,'lalpha',1)
   figure(1);   clf
   subplot(221);  box; hold on
   for f=1:length(Slopes.levels),
      if ~isempty(Slopes.up{f}) 
         plotedf(Slopes.up{f}); grid on
      end
      if ~isempty(Slopes.down{f}) 
         plotedf(-Slopes.down{f},'-.'); grid on 
      end
   end
   axis([0 10 0 1])
   title('\alpha = 0')
   subplot(222);  box; hold on
   for f=1:length(Slopes1.levels),
      if ~isempty(Slopes1.up{f}) 
         plotedf(Slopes1.up{f}); grid on
      end
      if ~isempty(Slopes1.down{f}) 
         plotedf(-Slopes1.down{f},'-.'); grid on 
      end
   end
   axis([0 10 0 1])
   title('\alpha = 1')
   subplot(223);  box; hold on
   plot(Slopes.meanwavex,Slopes.meanwaveup); 
   ax=axis;
   axis([Slopes.meanwavex(1) Slopes.meanwavex(end) ax(3) ax(4)]);
   plot(Slopes.meanwavex,Slopes.meanwavedown,'-.'); grid
   subplot(224);  box; hold on
   plot(Slopes1.meanwavex,Slopes1.meanwaveup);
   axis([Slopes.meanwavex(1) Slopes.meanwavex(end) ax(3) ax(4)]);
   plot(Slopes1.meanwavex,Slopes1.meanwavedown,'-.'); grid on
   pause
   
disp(' Section 1.3.2')
   % Figure: Comparison between empirical and theoretical slope CDF
   S = jonswap(1.5,[6 10]); S.h=20;
   opt = simoptset('Nt',2048*8,'dt',0.125,...
        'Nu',321,'du',0.125,'ffttype','ffttime');
   Nsim = 10;
   Slopes = spec2slopedat(S,Nsim,'time',[],opt);
   Slopes1 = spec2slopedat(S,Nsim,'time',[],opt,'lalpha',1);
   relativelevels = -1:2;
   y=0:0.01:10;
   [Fu,Fd]  = spec2timeslopecdf(S,y,relativelevels,opt);
   [Fu1,Fd1] = spec2timeslopecdf(S,y,relativelevels,opt,'lalpha',1);
   pause
   
   % Plot the results
   figure(1);   clf
   plotedf(Slopes1.up{4}); % Some very large values set the axis
   clf
   axis([0 10 0 1])
   hold on;    grid on
   plot(Fu1.x,Fu1.f{4},'r')
   plot(Fu1.x,Fu1.f{3},'r')
   plot(Fu1.x,Fu1.f{2},'r')
   plot(Fu1.x,Fu1.f{1},'r')
   plotedf(Slopes1.up{4});
   plotedf(Slopes1.up{3});
   plotedf(Slopes1.up{2});
   plotedf(Slopes1.up{1});  
   plot(Fd1.x,Fd1.f{1},'r-.')
   plot(Fd1.x,Fd1.f{2},'r-.')
   plot(Fd1.x,Fd1.f{3},'r-.')
   plot(Fd1.x,Fd1.f{4},'r-.')
   plotedf(-Slopes1.down{1},'-.');
   plotedf(-Slopes1.down{2},'-.');
   plotedf(-Slopes1.down{3},'-.');
   plotedf(-Slopes1.down{4},'-.');
   pause
   
disp(' Section 1.3.3')
   % Table: Asymmetry measures
   S = jonswap(1.5,[6 10]); S.h=20;
   opt = simoptset('Nt',2048*8,'dt',0.125,...
        'Nu',321,'du',0.125,'ffttype','ffttime');
   Nsim = 10;
   [Slope0,Steep0,~] = spec2steepdat(S,Nsim,'time',[],opt);
   [Slope05,Steep05,~] = spec2steepdat(S,Nsim,'time',[],opt,'lalpha',0.5);
   [Slope10,Steep10,~] = spec2steepdat(S,Nsim,'time',[],opt,'lalpha',1);
   [Slope15,Steep15,~] = spec2steepdat(S,Nsim,'time',[],opt,'lalpha',1.5);
   [Slope20,Steep20,~] = spec2steepdat(S,Nsim,'time',[],opt,'lalpha',2.0);

   disp('End of Chapter 1')
   
   TSTOP=cputime
   disp(['Total time = ' num2str(TSTOP-TSTART) ' sec']);
pause(pstate)

