function wafoinit
%WAFOINIT  setup all global variables of the WAFODEMO

% Since the calculations are in a very small extent dependent on eachother,
% most of the calculations are located in the wafofigXX.m files.
% Some minor calculations and all the initializations of variables are
% are done here.

% GLOBAL WAFOFIGNUM 
%
% wafoXXXX functions only recalculates the values of empty 
% variables ==> saves alot of computation time

% TODO % needs to complete the fatigue , cycles and Markov models figures

%revised pab april2004  
% fixed a bug: options to spec2sdat changed.  
%revised pab Feb2004
% By pab 28.01.2000

% Figure parameter:
%~~~~~~~~~~~~~~~~~~
global WAFOFIGNUM  
if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global pwdstr  wafomenulabels Jxn Nxr 

if isempty(pwdstr)
  pwdstr=pwd; % save path to where the demo was started
end
Nfigs=10;
if isempty(wafomenulabels) % automatic generation of menu labels
  cd(fullfile(waforoot,'papers','wafodemo'));  
  wafomenulabels=cell(Nfigs,1);
  for ix=1:Nfigs
    wafomenulabels{ix} = geth1line(['wafofig' num2str(ix)],1);   
  end
  cd(pwdstr)
end

if isempty(Nxr) % load if not already loaded 
  Nxr=load('gfaksr89.dat');
  %xn=xn(1:10000,:); % used for debugging 
end

if isempty(Jxn) % load if not already loaded 
  jfile='yura87.dat';
  Jxn=load(jfile);
  if strcmp(jfile,'yura87.dat'), %Selecting only a part of the data
    Jxn=Jxn(5000:55000,[1 3]);
    Jxn(:,2)=detrendma(Jxn(:,2),1500); % remove the trend
 %   Jxn=Jxn(5000:55000,:);
  end
  %xn=xn(1:10000,:); % used for debugging 
end



if 0 % this is not in paper but could be included
  if isempty(Jmap)
    disp('Loading map over the Japan Sea...')
    Jmap = load('japansea.dat');
  end
  if isempty(Nmap)
    disp('Loading map over the North Sea...')
    Nmap = load('northsea.dat');
  end
end


% Define Globals used (that the user also  might change)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global Fs                          % sampling frequency for torsethaugen
global St Hm0 Tp                   % torsethaugen output, input Fig 2
global xt Nsim Iseed               % spec2sdat output , input    Fig1
global Ste L                       % dat2spec output, input      Fig2
global fTt fTte u Np nit speed     % spec2thpdf output, input    Fig3  
global Tt rate                     % dat2wa output, input        Fig3
global ma mb sp                    % spreading input             Fig4

global fTcfAc NNp Nh Nnit Nspeed   % spec2thpdf output, input    Fig5  
global NVcf NHd Nrate              % dat2steep output, input     Fig5

global fTcfAcTc JNp Jh Jnit Jspeed % spec2thpdf output, input    Fig6  
global JTcf JAc Jind Jrate         % dat2steep output, input     Fig6

global kdeTt kernel hs L2          % kdebin input: kernel, smoothing  (Fig3)
global kdeVcfHd Nkernel Nhs NL2    % and transformation  parameter,   (Fig5)
global kdeTcfAcTc Jkernel Jhs JL2  % respectively                     (Fig6)

% Set torsethaugen default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(Hm0),     Hm0 = 6;                end
if isempty(Tp),       Tp = 8;                end
if isempty(Fs),       Fs = 0.95238095238095; end % sampling frequency
w = linspace(0,pi*Fs,257).';
% fig4
if isempty(sp),       sp = 15;               end
if isempty(ma),       ma = 5;                end
if isempty(mb),       mb = -2.5;             end


% Set spec2sdat and dat2spec default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(Nsim),   Nsim = 1000;             end % 
if isempty(Iseed), Iseed = 1000;             end % seed for simulation
if isempty(L),         L = 80;               end

% Set dat2wa and dat2steep default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(rate),   rate = 8;end %fig3
if isempty(Nrate), Nrate = 8;end %fig5
if isempty(Jrate), Jrate = 8;end %fig6

% Set spec2thpdf default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Fig3
if isempty(u),           u = 0;   end
if isempty(Np),         Np = 33;  end
if isempty(nit),       nit = 4;   end
if isempty(speed),   speed = 5;   end
% fig5
if isempty(NNp),        NNp = 33;  end
if isempty(Nnit),     Nnit = -1;  end
if isempty(Nspeed), Nspeed = 7;   end
if isempty(Nh),         Nh = linspace(0,8,31);;end
% fig6
if isempty(JNp),        JNp = 33;  end
if isempty(Jnit),     Jnit = 4;   end
if isempty(Jspeed), Jspeed = 5;   end
if isempty(Jh),         Jh = linspace(0,6,31);end


% Set kdebin default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Fig3
if isempty(kernel),     kernel = 'epan';  end
if isempty(L2),             L2 = .5;      end
%Fig5
if isempty(Nkernel),   Nkernel = 'epan';  end
if isempty(NL2),           NL2 = [.5 .5]; end
%Fig6
if isempty(Jkernel),   Jkernel = 'epan';  end
if isempty(JL2),           JL2 = [1 .5];  end % used
% if hs, Nhs Jhs is empty then kdebin uses the 
% default values calculated  from the data (this is default)



% Initialization:
%~~~~~~~~~~~~~~~~

if isempty(St),
  St=torsethaugen(w,[Hm0 ,Tp]);
end
if isempty(xt)
  xt=spec2sdat(St,[Nsim,1],[],Iseed);
end
if isempty(Ste)
  Ste=dat2spec(xt,L,[],[],0.95);  % estimated Spectrum
end






