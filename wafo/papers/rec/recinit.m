function recinit
%RECINIT  setup all global variables of the RECDEMO

% Since the calculations are very dependent on eachother,
% all the time consuming calculations are located here.
% Plotting and some minor calculations are left to the recfigXX.m
% functions 
%
% GLOBAL RECFIGNUM <=2 : globals defined and datasets are loaded
%                  >=3 : reconstruction of xn, extraction of V and H and
%                        distribution fitting to V and H
%
% Recinit only recalculates the values of empty 
% variables if RECFIGNUM>=3 ==> saves alot of computation time

% revised pab feb2005
% updated call to kdebin
% revised pab feb2004
% revised pab 18.04.2001
% - updated call to reconstruct due to changes in all transformation
%    estimation functions
% - updated call to dat2steep
% By pab 28.01.2000

% Figure parameter:
%~~~~~~~~~~~~~~~~~~
global RECFIGNUM  
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  clear global RECFIGNUM
  return
end

global pwdstr  recmenulabels recfilename xn map wind 
if isempty(pwdstr)
  pwdstr=pwd; % save path to where the demo was started
end

if isempty(recmenulabels) % automatic generation of menu labels
  cd(fullfile(waforoot,'papers','rec'));  
  recmenulabels=cell(12,1);
  for ix=1:13
    recmenulabels{ix} = geth1line(['recfig' num2str(ix)],1);   
  end
  cd(pwdstr)
end

if isempty(recfilename),
  recfilename='gfaks89.dat';
end

if isempty(xn) % load if not already loaded 
  xn=load(recfilename);
  %xn=xn(1:10000,:); % used for debugging 
end
dT=xn(2,1)-xn(1,1); % sampling

if isempty(map)
  disp('Loading map over the North Sea...')
  map = load('northsea.dat');
end
if isempty(wind)
  disp('Loading wind measurements from Statfjord A...')
  wind = load('sfa89.dat');
end

% Define Globals used
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
global inds zcrit dcrit ddcrit        % findoutliers output, input
global Nsim L csm1 csm2 param         % reconstruct input
global xr g g2 test0 tobs mu1o mu1oStd % reconstruct output                
global Sr m0 m2 Tm02 Hm0 Vrms Hrms    % Spectral density, moments and Sea state parameters
global V H rate                       % dat2steep output, input
global phat res noverlap              % dist2dfit output, input
global sphat CSMA CSMB                % dist2dsmfun2 output, input
global ft2 fkde                       % dist2dpdf2 and kdebin outputs
global kernel hs L2                   % kdebin input: kernel, smoothing parameter and transformation, respectively


% RECONSTRUCTION 
%~~~~~~~~~~~~~~~~

% Set findoutliers default values:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(zcrit),  zcrit  = 0.02;               end
% set dcrit corresponding to 5m/s and ddcrit corresponding to g (maximum acceleration of Stokes wave is g/2) 
if isempty(dcrit),  dcrit  = 5*dT;               end % dcrit=1.23 used in  article
if isempty(ddcrit), ddcrit = gravity(61)*dT^2;   end % ddcrit=1.5*1.23 used in  article

% Set reconstruct default values:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(Nsim),   Nsim   = 6;                  end 
if isempty(csm1),   csm1   = 0.95;               end
if isempty(csm2),   csm2   = 0.05;               end
if isempty(param),  param  = [-5 5 501];        end


%  

% DISTRIBUTION FITTING
%~~~~~~~~~~~~~~~~~~~~~

% set dat2steep default value
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(rate),   rate   = 8;                  end % interpolation rate of data set
                                                     %  before extraction of parameters
%
%Set dist2dfit default values:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(res),        res    = 0.2;  end
if isempty(noverlap); noverlap = 0;    end

%Set dist2dsmfun2 default values:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(CSMA),         CSMA = 0.95; end
if isempty(CSMB),         CSMB = 0.95; end
						     
% Set kdebin default values: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
if isempty(kernel),     kernel = 'epan';  end
if isempty(L2),             L2 = [.5 .5]; end
% if hs and L2 is empty then kdebin uses the 
% default values calculated  from the data (v,h)						     



% RECONSTRUCTION 
%~~~~~~~~~~~~~~~~

if RECFIGNUM>2,
  if isempty(inds)
    disp('Finding spurious points')
    inds = findoutliers(xn,zcrit,dcrit,ddcrit);
  end
  if isempty(xr)
    disp('Reconstruction ...')
    disp('This may take a while (20 min -> 1 hour depending on machine, inputs etc...)')
    [xr,g,g2,test0,tobs,mu1o,mu1oStd]=reconstruct(xn,inds,Nsim,L,[],...
	troptset('csm',csm1,'gsm',csm2,'param',param));
  
    % Estimate spectral density and moments
    Sr=dat2spec(xr,L);
    m0=trapz(Sr.w,Sr.S); % or alternatively spec2mom
    m2=trapz(Sr.w,Sr.S);
    Hm0=4*sqrt(m0);
    Tm02=2*pi*sqrt(m0/m2);
    Vrms=2*Hm0/Tm02;
    Hrms=Hm0/sqrt(2);
  end
end


% DISTRIBUTION FITTING:
%~~~~~~~~~~~~~~~~~~~~~~


if RECFIGNUM>6,				 
  if isempty(V) || isempty(H)
    % Not extracting the missing points   
    Ns=find(isnan(xn(:,2)), 1 ) ;
    if ~isempty(Ns), [V , H]=dat2steep(xr(1:Ns-1,:),rate,0);end
    Ne=find(isnan(xn(:,2)), 1, 'last' ); 
    if isempty(Ne),Ne=0;end
    [V2 , H2]=dat2steep(xr(Ne+1:end,:),rate,0);
    V=[V(:);V2(:)];H=[H(:);H2(:)];
  end
  if isempty(phat)
    disp('Fitting a 2D distribution to the data')
    phat=fitmargcnd2d(V/Vrms,H/Hrms,{'weibull','rayleigh'},[res,noverlap]);
  end
  if isempty(sphat)
    sphat=margcnd2dsmfun2(phat,linspace(0,3.5),[CSMA,CSMB ],[1 1 ]); % smooth the estimated parameter values
  end

  % Calculate the theoretical density of fitted model
  if isempty(ft2)				    
    ft2=pdfmargcnd2d(linspace(0,4),linspace(0,4),sphat);
    ft2.labx{1}='v';
    ft2.labx{2}='h';
  end
end


if RECFIGNUM>8, % kernel density estimation
  if isempty(fkde)
    disp('Estimating kernel density estimate')
    kopt = kdeoptset('kernel',kernel,'L2',L2,'inc',128);
    fkde = kdebin([V/Vrms,H/Hrms],kopt); % kdebin is much
                                                        % faster than kde
    k=find(fkde.f>30);
    if any(k)
      disp('Removing the spurios spikes')
      fkde.f(k)=0;
    end
    if 1,
      r = evalpdf(fkde,V/Vrms,H/Hrms,'linear');
      fkde.cl = qlevels2(r,fkde.pl); % calculate the levels which encloses fkde.pl
      % percent of the data (v,h)
    end
    fkde.labx{1}='v';
    fkde.labx{2}='h';
  end
end