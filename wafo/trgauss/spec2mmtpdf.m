function f = spec2mmtpdf(spec,utc,def,paramt,paramu,varargin) 
%SPEC2MMTPDF Joint density of Maximum, minimum and period.
%          
%  CALL:  f   = spec2mmtpdf(S,u,def,paramt,paramu,options); 
% 
%    f    = pdf (density structure) of crests (trough) heights 
%    S    = spectral density structure 
%    u    = reference level (default the most frequently crossed level).
%   def   = string defining density returned
%           'Mm'    : maximum and the following minimum. (M,m) (default)
%           'rfc'   : maximum and the rainflow minimum height.
%           'AcAt'  : (crest,trough) heights. 
%           'vMm'   : level v separated Maximum and minimum   (M,m)_v 
%           'MmTMm' : maximum, minimum and period between (M,m,TMm)
%           'vMmTMm': level v separated Maximum, minimum and period
%                     between (M,m,TMm)_v
%           'MmTMd' : level v separated Maximum, minimum and the period
%                     from Max to level v-down-crossing (M,m,TMd)_v.
%           'MmTdm' : level v separated Maximum, minimum and the period from 
%                     level v-down-crossing to min. (M,m,Tdm)_v
%           NB! All 'T' above can be replaced by 'L' to get wave length
%           instead.   
% paramt  = [0 tn Nt] defines discretization of half period: tn is the
%           longest period considered while Nt is the number of points,
%           i.e. (Nt-1)/tn is the sampling frequnecy. paramt= [0 10 51]
%           implies that the halfperiods are considered at 51 linearly
%           spaced points in the interval [0,10], i.e. sampling frequency
%           is 5 Hz.  
% paramu  = [u v N] defines discretization of maxima and minima ranges: 
%           u is the lowest minimum considered, v the highest maximum and N 
%           is the number of levels (u,v) included. 
% options = rind-options structure containing optional parameters
%           controlling the performance of the integration. 
%           See rindoptset for details.
%   []    = default values are used. 
%
% SPEC2MMTPDF calculates densities of wave characteristics in a
% stationary Gaussian transform process X(t) where 
%  Y(t) = g(X(t)) (Y zero-mean Gaussian with spectrum given in input spec). 
%  The tr. g can be estimated using lc2tr, dat2tr, hermitetr or ochitr.
% 
% Examples:% The joint density of zero separated Max2min cycles in time (a);
%          % in space (b); AcAt in time for nonlinear sea model (c): 
%
%    Hm0=7;Tp=11;
%    S = jonswap(4*pi/Tp,[Hm0 Tp]); 
%    Sk = spec2spec(S,'k1d');     
%    L0 = spec2mom(S,1); 
%    paramu = [sqrt(L0)*[-4 4] 41]; 
%    ft = spec2mmtpdf(S,0,'vmm',[],paramu); pdfplot(ft)         % a)
%    fs = spec2mmtpdf(Sk,0,'vmm');  figure, pdfplot(fs)         % b)
%    [sk, ku, me]=spec2skew(S); 
%    g = hermitetr([],[sqrt(L0) sk ku me]);
%    Snorm=S; Snorm.S=S.S/L0; Snorm.tr=g;
%    ftg=spec2mmtpdf(Snorm,0,'AcAt',[],paramu); pdfplot(ftg)    % c)
%
% See also  rindoptset, dat2tr, datastructures, wavedef, perioddef
 
% References
% Podgorski et al. (2000)
% "Exact distributions for apparent waves in irregular seas"
% Ocean Engineering,  Vol 27, no 1, pp979-1016.
% 
% P. A. Brodtkorb (2004), 
% Numerical evaluation of multinormal expectations
% In Lund university report series
% and in the Dr.Ing thesis: 
% The probability of Occurrence of dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
%
% Per A. Brodtkorb (2006)
% "Evaluating Nearly Singular Multinormal Expectations with Application to
% Wave Distributions",
% Methodology And Computing In Applied Probability, Volume 8, Number 1, pp. 65-91(27) 


% Tested on : matlab 5.3 
% History: by I. Rychlik 01.10.1998 with name minmax.m, new name: spec2cmat 
% bounds by I.R. 02.01.2000. 
% Revised by pab 09.05.2000 
%  - added level u separated Max and min + period distributions 
%  - new name: spec2cmat -> spec2mmtpdf 
%  - also Made call with directional spectrum possible.
% revised by IR removing error in transformation 29 VI 2000
% revised by IR  changed default value for paramt and some ch. of
% help 2 VII 2000
% revised by IR, addopted for Matlab 6, 31 III 2001 
% revised pab 10 April 2003
%  -added matlab version, i.e., a call to sp2mmtpdf  
%revised pab 22Nov 2003
% removed  nit and speed from input and replaced it 
%  with a rindoptions structure  
%    nit  =  0,...,9. Dimension of numerical integration  (default 2).
%           -1,-2,... different important sampling type integrations.
%   speed = defines accuracy of calculations, by choosing different 
%               parameters, possible values: 1,2...,9 (9 fastest, default
%               4).
  
defaultSpeed   = 4;
defaultoptions = rindoptset('speed',defaultSpeed,'nit',2,'method',0);
defaultoptions.speed = [];
  
if ((nargin==1) && (nargout <= 1) &&  isequal(spec,'defaults')),
  f = defaultoptions;
  return
end 
%error(nargchk(1,7,nargin))
narginchk(1,7)
startTime = clock; 

ftype = freqtype(spec);
if nargin<3||isempty(def) 
  defnr=0; 
  if strcmp(ftype,'k')
    def='L'; % distributions in space are default   
  else
    def='T'; % distributions in time are default   
  end
else
  switch lower(def)
    case 'acat',
      defnr =-2; % (Ac,At)
    case 'rfc',
      defnr =-1; % (M,m_rfc)
    case 'mm',
      defnr = 0; % max2min. (M,m) (default)
    case {'mmtmm','mmlmm'},
      defnr = 1; % max2min and period inbetween (M,m,TMm)
    case {'vmm'},
      defnr = 2; % level v separated Max2min   (Mv,mv)
    case {'vmmtmm',  'vmmlmm'},
      defnr = 3; % level v separated Max2min and period inbetween (Mu,mu,TMm)
    case {'mmtmd','vmmtmd','mmlmd','vmmlmd'},
      defnr = 4; % level v separated Max2min and period from Max to level v-down-crossing.
    case {'mmtdm','vmmtdm','mmldm','vmmldm'},
      defnr = 5; % level v separated Max2min and period from level v-down-crossing to min.
    otherwise, error('Unknown def')
  end
  if defnr>=3||defnr==1,
    def = upper(def(end-2)); % Store the kind of distribution that is
    % wanted: Period or wavelength
  elseif strcmp(ftype,'k')
    def='L'; % distributions in space are default
  else
    def='T'; % distributions in time are default
  end
end


switch upper(def(1))
  case {'L'},  spec = spec2spec(spec,'k1d') ;  ptxt='space';
  case {'T'},  spec = spec2spec(spec,'freq');  ptxt='time';
  otherwise,  error('Unknown def')
end

[S, xl4, L0, L2, L4] = specnorm(spec);
A    = sqrt(L0/L2); 

if isfield(spec,'tr') 
   g=spec.tr; 
else 
  g=[];
end 
if isempty(g)
  g  = [sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)'];
end

if nargin<6
  options = defaultoptions;
else
  options = rindoptset(defaultoptions,varargin{:});
end



if nargin<2||isempty(utc)
  utc_d = gaus2dat([0 0],g); % most frequent crossed level 
  utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0. utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ' num2str(u)])
	     

if nargin<4||isempty(paramt) 
  paramt = [0 5*pi*sqrt(L2/L4), 41];  
end 
if nargin<5||isempty(paramu)  
  paramu=[-5*sqrt(L0) 5*sqrt(L0) 41]; 
end 

t0     = paramt(1); 
tn     = paramt(2); 
Ntime  = paramt(3);

t      = levels([0 tn/A Ntime]); % normalized times 
Nstart = 1 + round(t0/tn*(Ntime-1)); % the starting point to evaluate 


Nx = paramu(3);
if (defnr>1)
  paramu(1) = max(0,paramu(1));
  if (paramu(2)<0), 
    error('Discretization levels must be larger than zero'), 
  end  
end
  
%Transform amplitudes to Gaussian levels:    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h   = levels(paramu); 
h   = reshape(h,Nx,1); 
der = ones(Nx,1);

if defnr>1 % level v separated Max2min densities
  der1 = der;
  hg   = tranproc([utc+h der],g);
  der  = abs(hg(:,2)); 
  hg   = hg(:,1); % Gaussian levels above u
  hg1  = tranproc([utc-h der1],g);
  der1 = abs(hg1(:,2)); 
  hg   = [hg;hg1(:,1)]; % Gaussian levels below u
else   % Max2min densities

  hg  = tranproc([h der],g); 
  der = abs(hg(:,2)); 
  der1= der;
  hg  = hg(:,1); % Gaussian level 
end

dt      = t(2)-t(1) ;

% Calculating covariances
%~~~~~~~~~~~~~~~~~~~~~~~~~
nr = 4; % number of derivatives
R = spec2cov2(S,nr,Ntime-1,dt);
 
%NB!!! the spec2XXpdf.exe programmes is very sensitive to how you interpolate  
%      the covariances, especially where the process is very dependent 
%      and the covariance matrix is nearly singular. (i.e. for small t 
%      and high levels of u if Tc and low levels of u if Tt) 
%     The best is to interpolate the spectrum linearly so that S.S>=0
%     This makes sure that the covariance matrix is positive
%     semi-definitt, since the circulant spectrum are the eigenvalues of
%     the circulant covariance matrix.


callFortran = 0; %options.method<0;  
if callFortran, % call fortran
   ftmp = cov2mmtpdfexe(R,dt,u,defnr,Nstart,hg,options);
   err = repmat(nan,size(ftmp));
else
  [ftmp,err,terr,options] = cov2mmtpdf(R,dt,u,defnr,Nstart,hg,options);   
end

f       = createpdf; 
if isfield(S,'note')
  f.note  = S.note;
end
f.options = options;
%f.u     = utc;

f.elapsedTime = etime(clock,startTime);

if Nx>2
  f.labx{1}='Max [m]';
  f.x{1}=h; 
  switch defnr 
    case -2, f.title = sprintf('Joint density of (Ac,At) in %s', ptxt);
    case -1, f.title = sprintf('Joint density of (M,m_{rfc}) in %s', ptxt);
    case  0, f.title = sprintf('Joint density of (M,m) in %s', ptxt);
    case  1,
      f.title = sprintf('Joint density of (M,m,%sMm) in %s', def(1), ptxt);
    case  2,
      f.title = sprintf('Joint density of (M,m)_{v=%2.5g} in %s',utc, ptxt);
    case  3,
      f.title = sprintf('Joint density of (M,m,%sMm)_{v=%2.5g} in %s',...
        def(1),utc, ptxt);
    case  4,
      f.title = sprintf('Joint density of (M,m,%sMd)_{v=%2.5g} in %s',...
        def(1),utc, ptxt);
    case  5,
      f.title = sprintf('Joint density of (M,m,%sdm)_{v=%2.5g} in %s',...
        def(1),utc, ptxt);
    otherwise, error('Unknown def')
  end
else
  f.note = [f.note 'Density is not scaled to unity'];
  switch defnr
    case {-2,-1,0,1},
      f.title = sprintf('Density of (%sMm, M = %2.5g, m = %2.5g)',...
        def(1),h(2),h(1));
    case {2,3},
      f.title = sprintf('Density of (%sMm, M = %2.5g, m = %2.5g)_{v=%2.5g}',...
        def(1),h(2),-h(2),utc);
    case 4,
      f.title = sprintf('Density of (%sMd, %sMm, M = %2.5g, m = %2.5g)_{v=%2.5g}',...
        def(1),def(1),h(2),-h(2),utc);
    case 5,
      f.title = sprintf('Density of (%sdm, %sMm, M = %2.5g, m = %2.5g)_{v=%2.5g}',...
        def(1),def(1),h(2),-h(2),utc);
    otherwise, error('Unknown def')
  end
end

if Nx>2  % amplitude distributions wanted
  f.x{2}    = h;
  f.labx{2} = 'min [m]';
  if defnr>1||defnr==-2 ,f.u  = utc;end   % save level u
  
  if defnr>2 || defnr==1
    ftmp = reshape(ftmp,Nx,Nx,Ntime);
    err  = reshape(err,Nx,Nx,Ntime);
    der0 = der1(:,ones(1,Nx)).*(der(:,ones(1,Nx)).');
    for ix =1:Ntime
      ftmp(:,:,ix) = ftmp(:,:,ix).*der0;% dh*dh;
      err(:,:,ix)  = err(:,:,ix).*der0;
    end
     
    ftmp   = ftmp/A;
    err    = err/A;
    f.x{3} = t(:)*A;
    if strcmpi(def(1),'t')
      f.labx{3} = 'period [sec]';
    else
      f.labx{3} = 'wave length [m]';
    end
  else
    ftmp  = reshape(ftmp,Nx,Nx);
    err   = reshape(err,Nx,Nx);
    ftmp  = ftmp.*der(:,ones(1,Nx)).*(der(:,ones(1,Nx)).');
    err   = err.*der(:,ones(1,Nx)).*(der(:,ones(1,Nx)).');
    if (defnr==-1)
      ftmp0 = fliplr(mctp2rfc(fliplr(ftmp)));
      err  = abs(ftmp0-fliplr(mctp2rfc(fliplr(ftmp+err))));
      ftmp = ftmp0;
    elseif (defnr==-2)
      ftmp0=fliplr(mctp2tc(fliplr(ftmp),utc,paramu))*sqrt(L4*L0)/L2;
      err =abs(ftmp0-fliplr(mctp2tc(fliplr(ftmp+err),utc,paramu))*sqrt(L4*L0)/L2);
      index1=find(f.x{1}>0);
      index2=find(f.x{2}<0);
      ftmp=flipud(ftmp0(index2,index1));
      err =flipud(err(index2,index1));
      f.x{1} = f.x{1}(index1);
      f.x{2} = abs(flipud(f.x{2}(index2)));
    end
  end
  f.f = ftmp;
  f.err = err;
else % Only time or wave length distributions wanted
  f.f = ftmp/A;
  f.err = err/A;
  f.x{1}=A*t';
  if strcmpi(def(1),'t')
    f.labx{1} = 'period [sec]';
  else
    f.labx{1} = 'wave length [m]';
  end
  if defnr>3,
    f.f   = reshape(f.f,[Ntime, Ntime]);
    f.err = reshape(f.err,[Ntime, Ntime]);
    f.x{2}= A*t';
    if strcmpi(def(1),'t')
      f.labx{2} = 'period [sec]';
    else
      f.labx{2} = 'wave length [m]';
    end
  end 
end
 

try
 [f.cl,f.pl]=qlevels(f.f,[10 30 50 70 90 95 99 99.9],f.x{1},f.x{2});
catch
  warning('WAFO:SPEC2MMTPDF','Singularity likely in pdf')
end
%pdfplot(f) 

%Test of spec2mmtpdf
% cd  f:\matlab\matlab\wafo\source\sp2thpdfalan
% addpath f:\matlab\matlab\wafo ,initwafo, addpath f:\matlab\matlab\graphutil
% Hm0=7;Tp=11; S = jonswap(4*pi/Tp,[Hm0 Tp]);
% ft = spec2mmtpdf(S,0,'vMmTMm',[0.3,.4,11],[0 .00005 2]);

return % main

function dens  = cov2mmtpdfexe(R,dt,u,defnr,Nstart,hg,options)
% Write parameters to file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Nx    = max(1,length(hg));
  if (defnr>1) 
    Nx = Nx/2; %level v separated max2min densities wanted     
  end
  Ntime = size(R,1);
  
  filenames = {'h.in','reflev.in'};
  cleanup(filenames{:}) 
  
  fid = fopen('h.in','wt');
  fprintf(fid,'%12.10f\n',hg);
  fclose(fid);
  
  %XSPLT = options.xsplit;
  nit   = options.nit;
  speed = options.speed;
  seed  = options.seed;
  SCIS  = abs(options.method); % method<=0
   
  disp('writing data') 
  fid=fopen('reflev.in','wt'); 
  fprintf(fid,'%2.0f \n',Ntime); 
  fprintf(fid,'%2.0f \n',Nstart); 
  fprintf(fid,'%2.0f \n',nit); 
  fprintf(fid,'%2.0f \n',speed); 
  fprintf(fid,'%2.0f \n',SCIS); 
  fprintf(fid,'%2.0f \n',seed);  % select a random seed for rind  
  fprintf(fid,'%2.0f \n',Nx); 
  fprintf(fid,'%12.10E \n',dt);
  fprintf(fid,'%12.10E \n',u);
  fprintf(fid,'%2.0f \n',defnr); % def
  fclose(fid);  
  
  filenames2 = writecov(R);
  
  disp('   Starting Fortran executable.') 

  dos([ wafoexepath 'cov2mmtpdf.exe']); %compiled cov2mmtpdf.f with rind70.f

  dens = load('dens.out');
  
  cleanup(filenames{:},filenames2{:})
  
  %% Clean up
  
  return
  
  
  function [pdf,err,terr options] = cov2mmtpdf(R,dt,u,def,Nstart,hg,options)
%COV2MMTPDF Joint density of Maximum, minimum and period.
%
% CALL  [pdf, err, options] = cov2mmtpdf(R,dt,u,def,Nstart,hg,options)
%
% pdf     = calculated pdf size Nx x Ntime
% err     = error estimate   
% terr    = truncation error
% options = requested and actual rindoptions used in integration.
% R       = [R0,R1,R2,R3,R4] column vectors with autocovariance and its
%          derivatives, i.e., Ri (i=1:4) are vectors with the 1'st to 4'th
%          derivatives of R0.  size Ntime x Nr+1
% dt      = time spacing between covariance samples, i.e., 
%           between R0(1),R0(2).
% u       = crossing level
% def     = integer defining pdf calculated:
%           0 : maximum and the following minimum. (M,m) (default)
%           1 : level v separated Maximum and minimum   (M,m)_v 
%           2 : maximum, minimum and period between (M,m,TMm)
%           3 : level v separated Maximum, minimum and period
%               between (M,m,TMm)_v
%           4 : level v separated Maximum, minimum and the period
%               from Max to level v-down-crossing (M,m,TMd)_v.
%           5 : level v separated Maximum, minimum and the period from 
%               level v-down-crossing to min. (M,m,Tdm)_v 
% Nstart  = index to where to start calculation, i.e., t0 = t(Nstart)
% hg      = vector of amplitudes length Nx or 0
% options = rind options structure defining the integration parameters 
%  
% COV2MMTPDF computes joint density of the  maximum and the following    
% minimum or level u separated maxima and minima + period/wavelength    
%
% For DEF = 0,1 : (Maxima, Minima and period/wavelength)
%         = 2,3 : (Level v separated Maxima and Minima and
%                  period/wavelength between them)
%  
% If Nx==1 then the conditional  density for  period/wavelength between
% Maxima and Minima given the Max and Min is returned
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% Y=  X'(t2)..X'(ts)..X'(tn-1)||X''(t1) X''(tn)|| X'(t1) X'(tn)  X(t1) X(tn) 
% = [       Xt                   Xd                    Xc            ]  
%                                                                       
% Nt = tn-2, Nd = 2, Nc = 4
%                                                                       
% Xt= contains Nt time points in the indicator function                 
% Xd=    "     Nd    derivatives in Jacobian
% Xc=    "     Nc    variables to condition on                          
%                                                                       
% There are 3 (NI=4) regions with constant barriers:                    
% (indI(1)=0);     for i\in (indI(1),indI(2)]    Y(i)<0.                
% (indI(2)=Nt)  ;  for i\in (indI(2)+1,indI(3)], Y(i)<0 (deriv. X''(t1)) 
% (indI(3)=Nt+1);  for i\in (indI(3)+1,indI(4)], Y(i)>0 (deriv. X''(tn))  
% 
%
% For DEF = 4,5 (Level v separated Maxima and Minima and
%                period/wavelength from Max to crossing)
%
% If Nx==1 then the conditional joint density for  period/wavelength
% between Maxima, Minima and Max to level v crossing given the Max and
% the min is returned
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Y=  X'(t2)..X'(ts)..X'(tn-1)||X''(t1) X''(tn) X'(ts)|| X'(t1) X'(tn)  X(t1) X(tn) X(ts)
% = [       Xt                      Xd                     Xc            ]
%                                                                         
% Nt = tn-2, Nd = 3, Nc = 5
%                                                                         
% Xt= contains Nt time points in the indicator function                      
% Xd=    "     Nd    derivatives                                             
% Xc=    "     Nc    variables to condition on                               
%                                                                            
% There are 4 (NI=5) regions with constant barriers:                        
% (indI(1)=0);     for i\in (indI(1),indI(2)]    Y(i)<0.                    
% (indI(2)=Nt)  ;  for i\in (indI(2)+1,indI(3)], Y(i)<0 (deriv. X''(t1))    
% (indI(3)=Nt+1);  for i\in (indI(3)+1,indI(4)], Y(i)>0 (deriv. X''(tn))
% (indI(4)=Nt+2);  for i\in (indI(4)+1,indI(5)], Y(i)<0 (deriv. X'(ts))    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  
% History
% Revised pab 22Nov2003
%  -new inputarguments  
%  -added err to output  
%Revised pab 03.04.2003
% -translated from f90 to matlab
%Revised pab 22.04.2000
% - added mean separated min/max + (Tdm, TMd) period distributions
% - added scis 

R0 = R(:,1);
R1 = R(:,2);
R2 = R(:,3);
R3 = R(:,4);
R4 = R(:,5);

Ntime = length(R0);

Nx0  = max(1,length(hg));
Nx1  = Nx0;
%Nx0 = Nx1; % just plain Mm
if (def>1) 
  Nx1 = Nx0/2;
 %  Nx0 = 2*Nx1;   % level v separated max2min densities wanted     
end
%disp(sprintf('def = %d',def))

% ***** The bound 'infinity' is set to 100*sigma *****
XdInf = 100.d0*sqrt(R4(1));
XtInf = 100.d0*sqrt(-R2(1));

Nc = 4;
NI = 4;
Nd = 2; 
%Mb = 1; 
%Nj = 0;

Nstart = max(2,Nstart); 

symmetry = 0;
isOdd = mod(Nx1,2);
if (def<=1) %THEN % just plain Mm 
   Nx = Nx1*(Nx1-1)/2;
   IJ = (Nx1+isOdd)/2;         
   if (hg(1)+hg(Nx1)==0 && (hg(IJ)==0 ||hg(IJ)+hg(IJ+1)==0) ) 
     symmetry=0;
     disp(' Integration region symmetric')
     % May save Nx1-isOdd integrations in each time step 
     % This is not implemented yet.
     %Nx = Nx1*(Nx1-1)/2-Nx1+isOdd
   end

   % CC = normalizing constant = 1/ expected number of zero-up-crossings of X' 
   %CC = 2*pi*sqrt(-R2(1)/R4(1)); 
   %  XcScale = log(CC)
   XcScale = log(2*pi*sqrt(-R2(1)/R4(1)));
else
   % level u separated Mm
   Nx = (Nx1-1)*(Nx1-1);
   if ( abs(u)<=eps && hg(1)+hg(Nx1+1)==0 && (hg(Nx1)+hg(2*Nx1)==0) )
     symmetry=0;
     disp(' Integration region symmetric')
     % Not implemented for DEF <= 3
     %IF (DEF.LE.3) Nx = (Nx1-1)*(Nx1-2)/2 
   end

   if (def>3) %THEN
      Nstart = max(Nstart,3);
      Nc = 5;
      NI = 5; 
      Nd = 3;
   end %ENDIF
   %CC= normalizing constant= 1/ expected number of u-up-crossings of X
   %CC = 2*pi*sqrt(-R0(1)/R2(1))*exp(0.5D0*u*u/R0(1)); 
   XcScale = log(2*pi*sqrt(-R0(1)/R2(1))) + 0.5D0*u*u/R0(1);
end %ENDIF

options.('xcscale') = XcScale;
opt0 = struct2cell(options);
%opt0 = opt0(1:10);
%seed = [];
%opt0 =  {SCIS,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS,seed,NIT1};



if (Nx>1)
   if ((def==0 || def==2)) %  (M,m) or (M,m)v distribution wanted
      asize = [Nx1,Nx1];
   else           
     % (M,m,TMm), (M,m,TMm)v  (M,m,TMd)v or (M,M,Tdm)v distributions wanted 
      asize = [Nx1,Nx1,Ntime];
   end 
elseif (def>3) %
   % Conditional distribution for (TMd,TMm)v or (Tdm,TMm)v given (M,m)  wanted 
   asize = [1,Ntime,Ntime];
else
  % Conditional distribution for  (TMm) or (TMm)v given (M,m) wanted
   asize = [1,1,Ntime];
end
   
% Initialization
%~~~~~~~~~~~~~~~~~
pdf = zeros(asize);
err = pdf;
terr = pdf;

BIG  = zeros(Ntime+Nc+1,Ntime+Nc+1);
ex   = zeros(1, Ntime + Nc + 1);          

%fxind = zeros(Nx,1);
xc    = zeros(Nc,Nx);
   
indI  = zeros(1,NI);
a_up  = zeros(1,NI-1);
a_lo  = zeros(1,NI-1); 


% INFIN =  INTEGER, array of integration limits flags:  size 1 x Nb   (in)
    %            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
    %            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
    %            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
    %            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
%INFIN = repmat(0,1,NI-1);
%INFIN(3)  = 1;
a_up(1,3) = +XdInf;
a_lo(1,1:2) = [-XtInf -XdInf];
if (def>3)
   a_lo(1,4) = -XtInf;   
end

IJ = 0 ;
if (def<=1) % THEN     % Max2min and period/wavelength
   for I=2:Nx1
      J = IJ+I-1;
      xc(3,IJ+1:J) =  hg(I);
      xc(4,IJ+1:J) =  hg(1:I-1).';
      IJ = J;
   end %do
else
   % Level u separated Max2min
   xc(Nc,:) = u;
   % Hg(1) = Hg(Nx1+1)= u => start do loop at I=2 since by definition we must have:  minimum<u-level<Maximum
   for i=2:Nx1                
      J = IJ+Nx1-1;
      xc(3,IJ+1:J) =  hg(i);              % Max > u
      xc(4,IJ+1:J) =  hg(Nx1+2:2*Nx1).';    % Min < u
      IJ = J;
   end %do   
end %IF
if (def <=3)
   h11 = fwaitbar(0,[],sprintf('Please wait ...(start at: %s)',datestr(now)));

   for Ntd = Nstart:Ntime 
      %Ntd=tn
      Ntdc = Ntd+Nc;
      Nt   = Ntd-Nd;
      indI(2) = Nt;
      indI(3) = Nt+1;
      indI(4) = Ntd;
      % positive wave period 
      BIG(1:Ntdc,1:Ntdc) = covinput(BIG(1:Ntdc,1:Ntdc),R0,R1,R2,R3,R4,Ntd,0);
      [fxind,err0,terr0] = rind(BIG(1:Ntdc,1:Ntdc),ex(1:Ntdc),a_lo,a_up,indI, ...
			    xc,Nt,opt0{:});
      
      
      %fxind  = CC*rind(BIG(1:Ntdc,1:Ntdc),ex(1:Ntdc),xc,Nt,NIT1,...
      %speed1,indI,a_lo,a_up);

      
      if (Nx<2) %THEN
         % Density of TMm given the Max and the Min. Note that the
         % density is not scaled to unity
         pdf(1,1,Ntd) = fxind(1);
         err(1,1,Ntd) = err0(1).^2;
         terr(1,1,Ntd) = terr0(1);
         %GOTO 100
      else
      
         IJ = 0;
         switch def
         case {-2,-1,0},  % joint density of (Ac,At),(M,m_rfc) or (M,m).
            for i = 2:Nx1 
               J = IJ+i-1;
               pdf(1:i-1,i,1) = pdf(1:i-1,i,1)+fxind(IJ+1:J).'*dt;%*CC;
               err(1:i-1,i,1) = err(1:i-1,i,1)+(err0(IJ+1:J).'*dt).^2;
               terr(1:i-1,i,1) = terr(1:i-1,i,1)+(terr0(IJ+1:J).'*dt);
               IJ = J;
            end %do 
         case  {1}  % joint density of (M,m,TMm)
            for  i = 2:Nx1
               J = IJ+i-1;
               pdf(1:i-1,i,Ntd) = fxind(IJ+1:J).';%*CC
               err(1:i-1,i,Ntd) = (err0(IJ+1:J).').^2;%*CC
               terr(1:i-1,i,Ntd) = (terr0(IJ+1:J).');%*CC
               IJ = J;
            end %do
         case {2},  % joint density of level v separated (M,m)v
            for  i = 2:Nx1
               J = IJ+Nx1-1;
               pdf(2:Nx1,i,1) = pdf(2:Nx1,i,1)+fxind(IJ+1:J).'*dt;%*CC;
               err(2:Nx1,i,1) = err(2:Nx1,i,1)+(err0(IJ+1:J).'*dt).^2;
               terr(2:Nx1,i,1) = terr(2:Nx1,i,1)+(terr0(IJ+1:J).'*dt);
               IJ = J;
            end %do 
         case {3} 
            % joint density of level v separated (M,m,TMm)v
            for  i = 2:Nx1 
               J = IJ+Nx1-1;
               pdf(2:Nx1,i,Ntd) = pdf(2:Nx1,i,Ntd)+fxind(IJ+1:J).';%*CC;
               err(2:Nx1,i,Ntd) = err(2:Nx1,i,Ntd)+(err0(IJ+1:J).').^2;
               terr(2:Nx1,i,Ntd) = terr(2:Nx1,i,Ntd)+(terr0(IJ+1:J).');
               IJ = J;
            end %do 
         end % SELECT   
      end %ENDIF
      waitTxt = sprintf('%s Ready: %d of %d',datestr(now),Ntd,Ntime);
      fwaitbar(Ntd/Ntime,h11,waitTxt);
      
   end %do
   close(h11);
   err = sqrt(err);
   %  goto 800
else
   %200 continue
   waitTxt = sprintf('Please wait ...(start at: %s)',datestr(now));
   h11 = fwaitbar(0,[],waitTxt);
   tnold= -1;
   for tn = Nstart:Ntime,
     Ntd  = tn+1;
     Ntdc = Ntd + Nc;
     Nt   = Ntd - Nd;
     indI(2) = Nt;
     indI(3) = Nt + 1;
     indI(4) = Nt + 2;
     indI(5) = Ntd;
     
     if (~symmetry) %IF (SYMMETRY) GOTO 300
        
       for ts = 2:tn-1,
         % positive wave period
         BIG(1:Ntdc,1:Ntdc) = covinput(BIG(1:Ntdc,1:Ntdc),R0,R1,R2,R3, ...
           R4,tn,ts,tnold);
	  
         [fxind,err0,terr0] = rind(BIG(1:Ntdc,1:Ntdc),ex(1:Ntdc), ...
           a_lo,a_up,indI,xc,Nt,opt0{:});
	  
         %tnold = tn;
         switch (def),
           case {3,4}
             if (Nx==1) %THEN
               % Joint density (TMd,TMm) given the Max and the min.
               % Note the density is not scaled to unity
               pdf(1,ts,tn) = fxind(1);%*CC
               err(1,ts,tn) = err0(1).^2;%*CC
               terr(1,ts,tn) = terr0(1);%*CC
             else
               % 4,  gives level u separated Max2min and wave period
               % from Max to the crossing of level u (M,m,TMd).
               IJ = 0;
               for  i = 2:Nx1,
                 J = IJ+Nx1-1;
                 pdf(2:Nx1,i,ts) = pdf(2:Nx1,i,ts)+ fxind(IJ+1:J).'*dt;
                 err(2:Nx1,i,ts) = err(2:Nx1,i,ts)+ (err0(IJ+1:J).'*dt).^2;
                 terr(2:Nx1,i,ts) = terr(2:Nx1,i,ts)+ (terr0(IJ+1:J).'*dt);
                 IJ = J;
               end %do
             end
           case  {5}
             if (Nx==1)
               % Joint density (Tdm,TMm) given the Max and the min.
               % Note the density is not scaled to unity
               pdf(1,tn-ts+1,tn) = fxind(1); %*CC
               err(1,tn-ts+1,tn) = err0(1).^2;
               terr(1,tn-ts+1,tn) = terr0(1);
             else
               % 5,  gives level u separated Max2min and wave period from
               % the crossing of level u to the min (M,m,Tdm).
                  
               IJ = 0;
               for  i = 2:Nx1
                 J = IJ+Nx1-1;
                 pdf(2:Nx1,i,tn-ts+1)=pdf(2:Nx1,i,tn-ts+1) + fxind(IJ+1:J).'*dt;%*CC;
                 err(2:Nx1,i,tn-ts+1)=err(2:Nx1,i,tn-ts+1) + (err0(IJ+1:J).'*dt).^2;
                 terr(2:Nx1,i,tn-ts+1)=terr(2:Nx1,i,tn-ts+1) + (terr0(IJ+1:J).'*dt);
                 IJ = J;
               end %do
             end
         end % SELECT
       end%         enddo
     else % exploit symmetry
       %300   Symmetry
       for ts = 2:floor(Ntd/2)
         % Using the symmetry since U = 0 and the transformation is
         % linear.
         % positive wave period
         BIG(1:Ntdc,1:Ntdc) = covinput(BIG(1:Ntdc,1:Ntdc),R0,R1,R2,R3, ...
           R4,tn,ts,tnold);
	    
         [fxind,err0] = rind(BIG(1:Ntdc,1:Ntdc),ex,a_lo,a_up,indI, ...
           xc,Nt,opt0{:});
	    
         %tnold = tn;
         if (Nx==1) % THEN
           % Joint density of (TMd,TMm),(Tdm,TMm) given the max and
           % the min.
           % Note that the density is not scaled to unity
           pdf(1,ts,tn) = fxind(1);%*CC;
           err(1,ts,tn) = err0(1).^2;
           err(1,ts,tn) = terr0(1);
           if (ts<tn-ts+1) %THEN
             pdf(1,tn-ts+1,tn) = fxind(1);%*CC;
             err(1,tn-ts+1,tn) = err0(1).^2;
             terr(1,tn-ts+1,tn) = terr0(1);
           end
           %GOTO 350
         else
           IJ = 0 ;
           switch (def)
             case {4}
               
               % 4,  gives level u separated Max2min and wave period from
               % Max to the crossing of level u (M,m,TMd).
               for  i = 2:Nx1
                 J = IJ+Nx1-1;
                 pdf(2:Nx1,i,ts) = pdf(2:Nx1,i,ts) + fxind(IJ+1:J)*dt;%*CC;
                 err(2:Nx1,i,ts) = err(2:Nx1,i,ts) + (err0(IJ+1:J)*dt).^2;
                 terr(2:Nx1,i,ts) = terr(2:Nx1,i,ts) + (terr0(IJ+1:J)*dt);
                 if (ts<tn-ts+1)
                   % exploiting the symmetry
                   pdf(i,2:Nx1,tn-ts+1) = pdf(i,2:Nx1,tn-ts+1)+fxind(IJ+1:J)*dt;%*CC
                   err(i,2:Nx1,tn-ts+1) = err(i,2:Nx1,tn-ts+1)+(err0(IJ+1:J)*dt).^2;
                   terr(i,2:Nx1,tn-ts+1) = terr(i,2:Nx1,tn-ts+1)+(terr0(IJ+1:J)*dt);
                 end
                 IJ = J;
               end %do
             case {5}
               % 5,   gives level u separated Max2min and wave period
               % from the crossing of level u to min (M,m,Tdm).
               for  i = 2:Nx1,
                 J = IJ+Nx1-1;
                 
                 pdf(2:Nx1,i,tn-ts+1)=pdf(2:Nx1,i,tn-ts+1)+fxind(IJ+1:J)*dt;
                 err(2:Nx1,i,tn-ts+1)=err(2:Nx1,i,tn-ts+1)+(err0(IJ+1:J)*dt).^2;
                 terr(2:Nx1,i,tn-ts+1)=terr(2:Nx1,i,tn-ts+1)+(terr0(IJ+1:J)*dt);
                 if (ts<tn-ts+1)
                   %*CC; % exploiting the symmetry
                   pdf(i,2:Nx1,ts) = pdf(i,2:Nx1,ts)+ fxind(IJ+1:J)*dt;
                   err(i,2:Nx1,ts) = err(i,2:Nx1,ts)+ (err0(IJ+1:J)*dt).^2;
                   terr(i,2:Nx1,ts) = terr(i,2:Nx1,ts)+ (terr0(IJ+1:J)*dt);
                 end %ENDIF
                 IJ = J;
               end %do
           end %END SELECT
         end
         %350
       end %do
     end
     waitTxt = sprintf('%s Ready: %d of %d',datestr(now),tn,Ntime);
     fwaitbar(tn/Ntime,h11,waitTxt);
      
     %400     print *,'Ready: ',tn,' of ',Ntime
   end %do
   close(h11);
   err = sqrt(err);
end % if

%Nx1,size(pdf) def  Ntime
if (Nx>1)% THEN
  IJ = 1;
  if (def>2 || def ==1)
    IJ = Ntime;
  end
  pdf = pdf(1:Nx1,1:Nx1,1:IJ);
  err = err(1:Nx1,1:Nx1,1:IJ);
  terr = terr(1:Nx1,1:Nx1,1:IJ);
else
  IJ = 1;
  if (def>3)
    IJ = Ntime;
  end
  pdf = squeeze(pdf(1,1:IJ,1:Ntime));
  err = squeeze(err(1,1:IJ,1:Ntime));
  terr = squeeze(terr(1,1:IJ,1:Ntime));
end
return
      
   
function BIG = covinput(BIG,R0,R1,R2,R3,R4,tn,ts,tnold)
%COVINPUT Sets up the covariance matrix  
%
% CALL BIG = covinput(BIG, R0,R1,R2,R3,R4,tn,ts)
%   
%  BIG = covariance matrix for X = [Xt,Xd,Xc] in spec2mmtpdf problems.
%     
% The order of the variables in the covariance matrix
% are organized as follows: 
% for  ts <= 1:
%    X'(t2)..X'(ts),...,X'(tn-1) X''(t1),X''(tn)  X'(t1),X'(tn),X(t1),X(tn) 
% = [          Xt               |      Xd       |          Xc             ]
%
% for ts > =2:
%    X'(t2)..X'(ts),...,X'(tn-1) X''(t1),X''(tn) X'(ts)  X'(t1),X'(tn),X(t1),X(tn) X(ts) 
% = [          Xt               |      Xd               |          Xc             ]
%
% where 
%
% Xt= time points in the indicator function
% Xd= derivatives
% Xc=variables to condition on
 
% Computations of all covariances follows simple rules: Cov(X(t),X(s)) = r(t,s),
% then  Cov(X'(t),X(s))=dr(t,s)/dt.  Now for stationary X(t) we have
% a function r(tau) such that Cov(X(t),X(s))=r(s-t) (or r(t-s) will give the same result).
%
% Consequently  Cov(X'(t),X(s))    = -r'(s-t)    = -sign(s-t)*r'(|s-t|)
%               Cov(X'(t),X'(s))   = -r''(s-t)   = -r''(|s-t|)
%               Cov(X''(t),X'(s))  =  r'''(s-t)  = sign(s-t)*r'''(|s-t|)
%               Cov(X''(t),X(s))   =  r''(s-t)   = r''(|s-t|)
% Cov(X''(t),X''(s)) =  r''''(s-t) = r''''(|s-t|)

if nargin<9||isempty(tnold)
  tnold = -1;
end

    
if (ts>1) % THEN
   shft = 1;
   N    = tn + 5 + shft;
   %Cov(Xt,Xc)
   %for
    i = 1:tn-2;
    %j = abs(i+1-ts)
    %BIG(i,N)  = -sign(R1(j+1),R1(j+1)*dble(ts-i-1)) %cov(X'(ti+1),X(ts)) 
    j = i+1-ts;
    tau = abs(j)+1;
    %BIG(i,N)  = abs(R1(tau)).*sign(R1(tau).*j.');
    BIG(i,N)  = R1(tau).*sign(j.');
   %end
   %Cov(Xc)
   BIG(N         ,N) =  R0(1);       % cov(X(ts),X(ts))
   BIG(tn+shft+3 ,N) =  R0(ts);      % cov(X(t1),X(ts))
   BIG(tn+shft+4 ,N) =  R0(tn-ts+1); % cov(X(tn),X(ts))
   BIG(tn+shft+1 ,N) = -R1(ts);      % cov(X'(t1),X(ts))
   BIG(tn+shft+2 ,N) =  R1(tn-ts+1); % cov(X'(tn),X(ts))
   %Cov(Xd,Xc)
   BIG(tn-1 ,N) =  R2(ts);      %cov(X''(t1),X(ts))     
   BIG(tn   ,N) =  R2(tn-ts+1); %cov(X''(tn),X(ts))
   
   %ADD a level u crossing  at ts
   
   %Cov(Xt,Xd)
   %for
      i = 1:tn-2;
      j = abs(i+1-ts);
      BIG(i,tn+shft)  = -R2(j+1); %cov(X'(ti+1),X'(ts)) 
   %end
   %Cov(Xd)  
   BIG(tn+shft,tn+shft) = -R2(1);  %cov(X'(ts),X'(ts))
   BIG(tn-1   ,tn+shft) =  R3(ts); %cov(X''(t1),X'(ts))
   BIG(tn     ,tn+shft) = -R3(tn-ts+1);  %cov(X''(tn),X'(ts))

   %Cov(Xd,Xc)
   BIG(tn+shft ,N       ) =  0.d0;        %cov(X'(ts),X(ts)) 
   BIG(tn+shft,tn+shft+1) = -R2(ts);      % cov(X'(ts),X'(t1))
   BIG(tn+shft,tn+shft+2) = -R2(tn-ts+1); % cov(X'(ts),X'(tn))
   BIG(tn+shft,tn+shft+3) =  R1(ts);      % cov(X'(ts),X(t1))
   BIG(tn+shft,tn+shft+4) = -R1(tn-ts+1); % cov(X'(ts),X(tn))
        
  if (tnold == tn)
     % A previous call to covinput with tn==tnold has been made
     % need only to update  row and column N and tn+1 of big:
     return
%      % make lower triangular part equal to upper and then return
%      for j=1:tn+shft
%         BIG(N,j)      = BIG(j,N)
%         BIG(tn+shft,j) = BIG(j,tn+shft)   
%      end
%      for j=tn+shft+1:N-1
%         BIG(N,j) = BIG(j,N)
%         BIG(j,tn+shft) = BIG(tn+shft,j)   
%      end
%      return
   end %if
   %tnold = tn;
else
  % N = tn+4;
   shft = 0;
end %if
          
if (tn>2)      
   %for
   i=1:tn-2;
   %cov(Xt)
   %   for j=i:tn-2
   %     BIG(i,j) = -R2(j-i+1)              % cov(X'(ti+1),X'(tj+1))
   %  end %do
   
   BIG(i,i) = toeplitz(-R2(i));  % cov(Xt) =   % cov(X'(ti+1),X'(tj+1))


   %cov(Xt,Xc)
   BIG(i      ,tn+shft+1) = -R2(i+1);         %cov(X'(ti+1),X'(t1))  
   BIG(tn-1-i ,tn+shft+2) = -R2(i+1);         %cov(X'(ti+1),X'(tn))  
   BIG(i      ,tn+shft+3) =  R1(i+1);         %cov(X'(ti+1),X(t1))  
   BIG(tn-1-i ,tn+shft+4) = -R1(i+1);         %cov(X'(ti+1),X(tn))  
  
  %Cov(Xt,Xd)
  BIG(i,tn-1)       = R3(i+1);          %cov(X'(ti+1),X''(t1))  
  BIG(tn-1-i,tn)    =-R3(i+1);          %cov(X'(ti+1),X''(tn)) 
  %end %do
end      
%cov(Xd)
BIG(tn-1  ,tn-1  ) = R4(1);
BIG(tn-1  ,tn    ) = R4(tn);     %cov(X''(t1),X''(tn))
BIG(tn    ,tn    ) = R4(1);

%cov(Xc)
BIG(tn+shft+3 ,tn+shft+3) = R0(1);        % cov(X(t1),X(t1))
BIG(tn+shft+3 ,tn+shft+4) = R0(tn);       % cov(X(t1),X(tn))
BIG(tn+shft+1 ,tn+shft+3) = 0.d0 ;        % cov(X(t1),X'(t1))
BIG(tn+shft+2 ,tn+shft+3) = R1(tn);       % cov(X(t1),X'(tn))
BIG(tn+shft+4 ,tn+shft+4) = R0(1) ;       % cov(X(tn),X(tn))
BIG(tn+shft+1 ,tn+shft+4) =-R1(tn);       % cov(X(tn),X'(t1))
BIG(tn+shft+2 ,tn+shft+4) = 0.d0;         % cov(X(tn),X'(tn)) 
BIG(tn+shft+1 ,tn+shft+1) =-R2(1);        % cov(X'(t1),X'(t1))
BIG(tn+shft+1 ,tn+shft+2) =-R2(tn);       % cov(X'(t1),X'(tn))
BIG(tn+shft+2 ,tn+shft+2) =-R2(1);        % cov(X'(tn),X'(tn))
%Xc=X(t1),X(tn),X'(t1),X'(tn) 
%Xd=X''(t1),X''(tn)
%cov(Xd,Xc)
BIG(tn-1  ,tn+shft+3) = R2(1);           %cov(X''(t1),X(t1))     
BIG(tn-1  ,tn+shft+4) = R2(tn);          %cov(X''(t1),X(tn))     
BIG(tn-1  ,tn+shft+1) = 0.d0  ;          %cov(X''(t1),X'(t1))     
BIG(tn-1  ,tn+shft+2) = R3(tn);          %cov(X''(t1),X'(tn))     
BIG(tn    ,tn+shft+3) = R2(tn);          %cov(X''(tn),X(t1))     
BIG(tn    ,tn+shft+4) = R2(1);           %cov(X''(tn),X(tn))     
BIG(tn    ,tn+shft+1) =-R3(tn);          %cov(X''(tn),X'(t1))     
BIG(tn    ,tn+shft+2) = 0.d0;            %cov(X''(tn),X'(tn))
      
     
% make lower triangular part equal to upper 
%for j=1:N-1
%   for i=j+1:N
%      BIG(i,j) = BIG(j,i)
%   end %do
%end %do
lp      = find(tril(ones(size(BIG)),-1)); % indices to lower triangular part
BIGT    = BIG.';
BIG(lp) = BIGT(lp);
return
% END  SUBROUTINE COV_INPUT
        









  