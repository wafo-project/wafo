function recfig(ix)
% RECFIG  callback implementing functions of RecDemo
global RECFIGNUM 
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
Nfigs=13;
global recmenulabels


switch ix
  case -1,  % change default settings
    RECFIGNUM=0;
    changesettingsmenu; % local function    
  case -2, % make all figures
    RECFIGNUM=Nfigs;
    recinit
    cfig=gcf;
    for iy=1:Nfigs,
      figure(cfig-1+iy)
      figname = ['recfig',num2str(iy)];
      %clc; home;
      help(figname);
      eval(figname);
      drawnow;
    end
    for iy=1:Nfigs,
      figure(cfig-1+iy)
      figname = ['recfig',num2str(iy)];
      help(figname);
    end
   case -3, % select a figure
     if isempty(recmenulabels), recinit; end
     iy = menu('Choose a figure',recmenulabels{:},'Cancel');
     if ~isempty(iy) & (1 <= iy & iy <=Nfigs) ,
       RECFIGNUM = iy;
       recinit
       figname = ['recfig',num2str(iy)];
       
       clc; home; help(figname);
       eval(figname);
       drawnow;
     end
   otherwise
     disp('Unknown Input argument to recfig') 
     disp(sprintf('num =%g',ix))
end
     
   
function changesettingsmenu
% Local function implementing 
%

global RECFIGNUM  recfilename xn xr ft2 fkde V H
%Reconstruction parameters
global zcrit dcrit ddcrit Nsim L csm1 csm2 
% distribution fitting parameters
global rate  res noverlap CSMA CSMB
% kernel density estimation parameters
global kernel hs L2

RECFIGNUM=0;
recinit % initialize global variables
        % and put them into the cellarray of default values 
def=cell(3,1);
def{1}={zcrit,dcrit,ddcrit,Nsim,L,csm1,csm2};
def{2}={rate,res,noverlap,CSMA,CSMB};
def{3}={kernel,hs,L2};
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

answer=inputdlg({'Enter name of data set to reconstruct:'},'Change data set',lineNo,{recfilename},'on');
Na=length(answer);
if Na>0
  if (~strcmp(deblank(answer{ix}),deblank(recfilename))),
    recfilename = answer{ix};
    xn=[];
    xr=[];V=[];H=[];phat=[];sphat=[]; % all must be calculated again
    ft2=[];fkde=[];inds=[];
    sphat=[]; ft2=[];           
    phat=[];V=[];H=[]; fkde=[];
  end
end  

title1='Change Reconstruction parameters (See findoutliers, reconstruct)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt1 ={ ...
'Enter critical distance between consecutive points. (zcrit):',...
'Enter critical distance of diff(xn(:,2)) (dcrit):',...
'Enter critical distance of diff(xn(:,2),2) (ddcrit):',...
'Enter maximum number of iterations of the reconstruction (Nsim):', ...
'Enter lag size of the window function (default the lag where  ACF<2 stdev, maximum 200):',...
'Enter a value for csm1 which defines smoothing  of log(crossing intensity):',...
'Enter a value for csm2 which defines smoothing  of transformation g:'};
answer=inputdlg(prompt1,title1,lineNo,strdef{1},'on');
Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1>0); 
  for ix=ind2(:).',
    answer{ix}=str2num(answer{ix});
    if (isempty(answer{ix})|answer{ix}==def{1}{ix}), ind1(ix)=0; else, def{1}{ix}=answer{ix}; end
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    xr=[];V=[];H=[];phat=[];sphat=[]; % some or all of def{1} is changed
    ft2=[];fkde=[];				      
    if (any(ind2<4)),inds=[];end      % some or all of zcrit,dcrit or ddcrit is changed
    [zcrit,dcrit,ddcrit,Nsim,L,csm1,csm2]=deal(def{1}{:});
  end
end



title2='Change distribution fitting parameters (See dat2steep dist2dfit dist2dsmfun2)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
prompt2 ={...
'Enter interpolation rate of the data before extracting data (rate):',...
'Enter a value for res  which defines the bin size used in the conditional fitting:',...
'Enter a value for noverlap which defines the number of groups overlapping in the conditional fitting:',...
'Enter a value for CSMA which defines smoothing of Weibull scale parameter:',...
'Enter a value for CSMB which defines smoothing of Weibull shape parameter:'};
answer=inputdlg(prompt2,title2,lineNo,strdef{2},'on');	
Na=length(answer);
if Na>0
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1>0);
  for ix=ind2(:).',
    answer{ix}=str2num(answer{ix});
    if (isempty(answer{ix})|answer{ix}==def{2}{ix}), ind1(ix)=0; else, def{2}{ix}=answer{ix}; end
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    sphat=[]; ft2=[];                % some or all of def{2} is changed
    if (any(ind2<2)), phat=[];V=[];H=[]; fkde=[];end % rate is changed
    if (any(ind2==2 | ind2==3)), phat=[];        end % res or  noverlap  changed
    [rate,res,noverlap,CSMA,CSMB]=deal(def{2}{:});
  end
end



title3='Change Kernel Density Estimation parameters (See kdebin)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prompt3={...
'Enter name for kernel:',...
'Enter a smoothing parameter vector for hs (default 1D optimal value using hns):',...
'Enter vector of transformation parameters (L2):'};
answer=inputdlg(prompt3,title3,lineNo,strdef{3},'on');	
Na=length(answer);
if Na>0,
  ind1=zeros(1,Na); % index to changed variables
  for ix=1:Na, ind1(ix)=~isempty(answer(ix));  end
  ind2=find(ind1>0);
  for ix=ind2(:).',
    if (ix~=1), % ix==1 ==> kernel string nothing to do
      answer{ix}=str2num(answer{ix});
      if (isempty(answer{ix})|answer{ix}==def{3}{ix}), ind1(ix)=0; else, def{3}{ix}=answer{ix}; end
    else
      if (strcmpi(answer{ix},def{3}{ix})), ind1(ix)=0; else, def{3}{ix}=answer{ix}; end
    end
  end
  ind2=find(ind1>0); 
  if any(ind2) % set variables  that must be calculated again to empty.
    fkde=[]; % some or all of def{3} is changed	
    [kernel,hs,L2]=deal(def{3}{:});
  end
end