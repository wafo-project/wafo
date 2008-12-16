function wafofig(ix)
% WAFOFIG  callback implementing functions of WAFODEMO
global WAFOFIGNUM 
if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  return
end
Nfigs=10; % 
global wafomenulabels


switch ix
  case -1,  % change default settings
    WAFOFIGNUM=0;
    changesettingsmenu; % local function
    
  case -2, % make all figures
    WAFOFIGNUM=Nfigs;
    wafoinit
    cfig=gcf;
    for iy=1:Nfigs,
      figure(cfig-1+iy)
      figname = ['wafofig',num2str(iy)];
       %clc; home;
      help(figname);
      eval(figname); 
      drawnow;
    end
    for iy=1:Nfigs,
      figure(cfig-1+iy)
      figname = ['wafofig',num2str(iy)];
      help(figname);
    end
   case -3, % select a figure
     if isempty(wafomenulabels)| WAFOFIGNUM==0, wafoinit; end
     iy = menu('Choose a figure',wafomenulabels{:},'Cancel')
     if ~isempty(iy) & (1 <= iy & iy <=Nfigs) ,
       WAFOFIGNUM = iy;
       
       figname = ['wafofig',num2str(iy)];
       
       clc; home; help(figname);
       eval(figname);
       clc; home; help(figname);
       drawnow;
     end
     
   otherwise
     disp('Unknown Input argument to wafofig') 
     disp(sprintf('num =%g',ix))
end
     
   
function changesettingsmenu
%
%

global WAFOFIGNUM 
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
global fTcfAcTc JNp Jh Jnit Jspeed      % spec2thpdf output, input    Fig6  
global JTcf JAc Jind Jrate         % dat2steep output, input     Fig6
global kdeTt kernel hs L2          % kdebin input: kernel, smoothing  (Fig3)
global kdeVcfHd Nkernel Nhs NL2    % and transformation  parameter,   (Fig5)
global kdeTcfAcTc Jkernel Jhs JL2  % respectively                     (Fig6)

WAFOFIGNUM=0;
wafoinit % initialize global variables
        % and put them into the cellarray of default values 
def=cell(5,1);
def{1}={Fs Hm0 Tp sp ma mb};
def{2}={Nsim Iseed L };
def{3}={rate kernel hs L2 u Np nit speed};
def{4}={Nrate Nkernel Nhs NL2 NNp Nnit Nspeed};
def{5}={Jrate Jkernel Jhs JL2 JNp Jnit Jspeed};

% Convert def to strdef containing only strings
strdef=def;
for ix=1:length(def)
  for iy=1:length(def{ix})
    if isnumeric(def{ix}{iy})
      strdef{ix}{iy}=num2str(def{ix}{iy});
    end
  end
end

lineNo=1;
title1='Change Spectral parameters fig1-4 (See torsethaugen, spreading)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt1 ={ ...
'Enter sampling frequency, Fs. (Hz):',...
'Enter significant wave height, Hm0:',...
'Enter peak period, Tp:',...
'Enter spreading parameter, sp:', ...
'Enter spreading shape parameter, ma):',...
'Enter spreading shape parameter, mb):'};
answer=inputdlg(prompt1,title1,lineNo,strdef{1},'on');

Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1); 
  for ix=ind2(:).',
    answer{ix}=str2num(answer{ix});
    if (isempty(answer{ix})|answer{ix}==def{1}{ix}), ind1(ix)=0; else, def{1}{ix}=answer{ix}; end
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    if any(ind2<=3), % fig1,2,3
      St=[];xt=[];Ste=[]; Tt=[]; 
      fTt=[]; fTte=[];	kdeTt=[];
    end 
    [Fs, Hm0, Tp, sp, ma, mb]=deal(def{1}{:});
  end
end


title2='Change Simulation parameters fig1-3 (See spec2sdat, dat2spec, dat2wa)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt2 ={ ...
'Enter length of simulated time series. (Nsim):',...
'Enter seed for random number generator (iseed):',...
'Enter maximum lag size of the window function (L).:'};
answer=inputdlg(prompt2,title2,lineNo,strdef{2},'on');

Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1); 
  for ix=ind2(:).',
    answer{ix}=str2num(answer{ix});
    if (isempty(answer{ix})|answer{ix}==def{2}{ix}), ind1(ix)=0; else, def{2}{ix}=answer{ix}; end
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    Ste=[];fTte=[];
    if any(ind2<=2), % fig1,2,3
      xt=[]; Tt=[]; 
      kdeTt=[];
    end
    [Nsim, Iseed, L]=deal(def{2}{:});
  end
end

title3='Change distribution fitting parameters fig3 (See dat2wa kdebin spec2thpdf)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt3 ={ ...
'Enter interpolation rate before extracting Tt from data. (rate):',...
'Enter name for kernel:',...
'Enter smoothing parameter hs (default 1D optimal value using hns):',...
'Enter transformation parameter (L2):',...
'Enter reference level defining Tt (u):',...
'Enter number of points for which the theoretical pdf is computed (Np<100):', ......
'Enter nit (-2 and 0,1,...,9):',....
'Enter speed (1,2,...,9):'};

answer=inputdlg(prompt3,title3,lineNo,strdef{3},'on');

Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1); 
  for ix=ind2(:).',
    if ix~=2
      answer{ix}=str2num(answer{ix});
      if (isempty(answer{ix})|answer{ix}==def{3}{ix}), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    else
      if (strcmpi(answer{ix},def{3}{ix})), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    end
    
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    if any(ind2<=5), % fig3
     
      kdeTt=[];
      if any(ind2==1), Tt=[];  end
    end
    if any(ind2>=5), % fig3
       fTt=[]; fTte=[];
    end
    [rate kernel hs L2 u Np nit speed]=deal(def{3}{:});
  end
end


title4='Change distribution fitting parameters fig5 (See dat2steep kdebin spec2thpdf)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt4 ={ ...
'Enter interpolation rate before extracting parameters from data. (rate):',...
'Enter name for kernel:',...
'Enter smoothing parameter vector hs (default 1D optimal value using hns):',...
'Enter transformation parameter vector (L2):',...
'Enter number of points for time axis for which theoretical pdf is computed (Np<100):',...
'Enter nit (-2 and 0,1,...,9):',...
'Enter speed (1,2,...,9):'};
answer=inputdlg(prompt4,title4,lineNo,strdef{4},'on');

Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1); 
  for ix=ind2(:).',
    if ix~=2
      answer{ix}=str2num(answer{ix});
      if (isempty(answer{ix})|answer{ix}==def{3}{ix}), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    else
      if (strcmpi(answer{ix},def{3}{ix})), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    end
    
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    if any(ind2<=4), % fig5
      kdeVcfHd=[];
      if any(ind2==1) ,NVcf=[],NHd=[];end
    end
    if any(ind2>=4), % fig5
       fTcfAc=[]; 
    end
    [Nrate Nkernel Nhs NL2 NNp Nnit Nspeed]=deal(def{4}{:});
  end
end


title5='Change distribution fitting parameters fig6 (See dat2steep kdebin spec2thpdf)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt5 ={ ...
'Enter interpolation rate before extracting parameters from data. (rate):',...
'Enter name for kernel:',...
'Enter smoothing parameter vector hs (default 1D optimal value using hns):',...
'Enter transformation parameter vector (L2):',...
'Enter number of points for time axis for which theoretical pdf is computed (Np<100):',...
'Enter nit (-2 and 0,1,...,9):',...
'Enter speed (1,2,...,9):'};
answer=inputdlg(prompt5,title5,lineNo,strdef{5},'on');

Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1); 
  for ix=ind2(:).',
    if ix~=2
      answer{ix}=str2num(answer{ix});
      if (isempty(answer{ix})|answer{ix}==def{3}{ix}), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    else
      if (strcmpi(answer{ix},def{3}{ix})), ind1(ix)=0; else,def{3}{ix}=answer{ix}; end
    end
    
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    if any(ind2<=4), % fig6
      kdeTcfAcTc=[];
      if any(ind2==1) ,JTcf=[],JAc=[];end
    end
    if any(ind2>=4), % fig6
       fTcfAcTc=[]; 
    end
    [Jrate Jkernel Jhs JL2 JNp Jnit Jspeed]=deal(def{5}{:});
  end
end