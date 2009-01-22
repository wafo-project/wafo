function f = spec2thpdf(spec,utc,def,paramt,h,options,plotflag)
%SPEC2THPDF Joint density of amplitude and period/wave-length characteristics
%         
%  CALL:  f   = spec2thpdf(S,u,def,paramt,h,options,plotflag);
%
%    f    = density structure of wave characteristics of half-waves
%    S    = spectral density structure
%    u    = reference level (default the most frequently crossed level).
%    def  = 'Tc',    gives crest period, Tc (default).
%           'Tt',    gives trough period, Tt.
%           'TcAc',  gives crest period and wave crest amplitude (Tc,Ac).
%           'TtAt',  gives trough period and wave trough amplitude (Tt,At).
%           'TcfAc', gives crest front period and wave crest amplitude (Tcf,Ac).
%           'TtbAt', gives trough back period and wave trough amplitude (Ttb,At).
%           'TAc',   gives minimum of crest front/back period and wave crest 
%                                                  amplitude (min(Tcf,Tcb),Ac).
%           'TAt',   gives minimum of trough front/back period and wave trough 
%                                                 amplitude (min(Ttf,Ttb),At). 
%            NB! All 'T.' above can be replaced by 'L.' to get wave
%            length instead.
% paramt  = [t0 tn Nt] where t0, tn and Nt is the first value, last value 
%               and the number of points, respectively, for which
%               the density will be computed. paramt= [5 5 51] implies
%               that the density is computed only for T=5 and
%               using 51 linearly spaced points in the interval [0,5].
%    h    = vector of  amplitudes (default levels([0 hmax 31]))  note  h >= 0
%           hmax evaluated from spec
% options = rind-options structure containing optional parameters
%           controlling the performance of the integration.
%           See rindoptset for details.
%plotflag = plots result if 1, no plot else (default 0)  
%      [] = default values are used.
%
% SPEC2THPDF calculates densities of wave characteristics of half waves 
% in a stationary Gaussian transform process X(t) where 
% Y(t) = g(X(t)) (Y zero-mean Gaussian with spectrum given in input spec). 
% The transformation, g, can be estimated or calculated using LC2TR,
% DAT2TR, HERMITETR or OCHITR.  
% 
% Note that due to symmetry the density of (Tcf,Ac) is equal to (Tcb,Ac).
% Similarly  the density of (Ttf,At) is equal to (Ttb,At).
%
% Example:  The density of (Tcf,Ac,Tc=5) is computed by
%
%    opt = rindoptset('speed',5,'nit',3,'method',0);
%  
%    S = jonswap;
%    f = spec2thpdf(S,[],'TcfAc',[5 5 51],[],opt,1);
%    
%    opt1 = rindoptset( 'speed',9,'nit',3,'method',1); 
%    paramt = [0 10,51] 
%    f = spec2thpdf(S,[],'Tc',paramt,[],opt1);
%    pdfplot(f)
%
% See also   rindoptset, dat2tr, datastructures, ampdef, perioddef

% Tested on : matlab 5.3,6.X
% History: by I. Rychlik 01.10.1998 with name wave_th1.m
% revised by Per A. Brodtkorb 19.09.1999
% revised by I.R. 30.09.1999, bugs removing.
% revised by pab 21.10.1999
% added option 'TAc', 'TAt'
% updated to new pdf structure 
% revised by es 000322.  Make call with directional spectrum possible.
%                        Introduced 'L..' options for def.  
%                        Added plotflag and made some changes in the help  
% revised pab 28.03.2000
%  - added NIT=-3:-11
%  - added XSPLT, rateLHD to reflev.in
%  - fixed 2 bugs: 1) contour levels only for 2D pdf's
%                  2) title string missed some braces
% revised pab 22.05.2000
%  - changed order of negative NIT's ie SCIS

%               -1 Integrate all by SADAPT for Ndim<9 and by KRBVRC otherwise 
%               -2 Integrate all by SADAPT by Genz (1992) (Fast)
%               -3 Integrate all by KRBVRC by Genz (1993) (Fast)
%               -4 Integrate all by KROBOV by Genz (1992) (Fast)
%               -5 Integrate all by RCRUDE by Genz (1992)
%               -6 Integrate all by RLHCRUDE using MLHD and center of cell 
%               -7 Integrate all by RLHCRUDE using  LHD and center of cell
%               -8 Integrate all by RLHCRUDE using MLHD and random point within the cell
%               -9 Integrate all by RLHCRUDE using LHD and random point within the cell
% revised by IR removing error in transformation 29 VI 2000
% revised by IR removing error in normalization when t0=tn defnr=3,-3 1 VII 2000
% revised by IR adapted to MATLAB6 + reshaping h, 8 II 2001
% revised by jr 01.03.28  changed 'h.in' to fid at line 181
% revised pab 03.03.03 
%  - added sp2thpdf.m  
% revised pab 19.11.2003  
% -added err to output
% -removed some code and replaced it with a call to spec2cov2  
% -removed nit and speed from input/output and replaced it 
%  with a rindoptions structure  
%    nit  =  0,...,9. Dimension of numerical integration  (default 2).
%           -1,-2,... different important sampling type integrations.
%   speed = defines accuracy of calculations, by choosing different 
%               parameters, possible values: 1,2...,9 (9 fastest, default 4). 
% Revised pab Aug2007
% -replaced some code with call to writecov

  
defaultSpeed   = 9;
defaultoptions = rindoptset('speed',defaultSpeed);
%defaultoptions.nit = 2;
defaultoptions.speed = [];
defaultoptions.coveps = 0.01;
  
if ((nargin==1) && (nargout <= 1) &&  isequal(spec,'defaults')),
  f = defaultoptions;
  return
end 
error(nargchk(1,7,nargin))
startTime = clock;
if nargin<3||isempty(def)
  def='tc';
end
if strcmpi('l',def(1))
  spec = spec2spec(spec,'k1d');
elseif strcmpi('t',def(1))
  spec = spec2spec(spec,'freq');
else
  error('Unknown def')
end
switch lower(def(2:end))
 case  'c',               defnr = 1;
 case  't',               defnr =-1;
 case  'cac',             defnr = 2;
 case  'tat',             defnr =-2;
 case  {'cfac','cbac'} , defnr = 3;
 case  {'tbat', 'tfat'}, defnr =-3;
 case  {'ac'} ,          defnr = 4;
 case  {'at'},           defnr =-4;
 otherwise, error('Unknown def')
end


%[S, xl4, L0, L2, L4, L1] = specnorm(spec);

[S, xl4, L0, L2] = specnorm(spec);
A   = sqrt(L0/L2);

if isfield(spec,'tr')
   g = spec.tr;
else
   g = [];
end
if isempty(g)
  g  = [sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)'];
end
if nargin<2||isempty(utc)
  utc_d = gaus2dat([0, 0],g); % most frequently crossed level 
  utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0., utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ', num2str(u)])

if nargin<6||isempty(options)
  options = defaultoptions;
else
  options = rindoptset(defaultoptions,options);
end

if nargin<7||isempty(plotflag)
  plotflag=0;
end			    

if nargin<4||isempty(paramt)
  %z2 = u^2/2;
  z  = -sign(defnr)*u/sqrt(2);
  expectedMaxPeriod = 1.5*ceil(2*pi*A*exp(z)*(0.5+erf(z)/2)); 
  paramt = [0, expectedMaxPeriod, 51];
end

t0     = paramt(1);
tn     = paramt(2);
Ntime  = paramt(3);
t    = levels([0, tn/A, Ntime]); % normalized times
Nstart = 1 + round(t0/tn*(Ntime-1)); % index to starting point to evaluate

if abs(defnr)<2 
  nr = 2;
  Nx = 1;
  hg = [];
else
  nr = 4;
  if nargin<5||isempty(h) 
    Nx = 31;
    px = gaus2dat([0., u;1, 4.*sign(defnr)],g); 
    px = abs(px(2,2)-px(1,2));
    h  = levels([0, px, Nx]).';
  else
    h  = sort(abs(h(:))); % make sure values are positive and increasing
    Nx = length(h);
    h  = reshape(h,Nx,1);
  end
  %Transform amplitudes to Gaussian levels:   
  der = ones(Nx,1); % dh/dh=1
  hg  = tranproc([utc+sign(defnr)*h, der],g);
  der = abs(hg(:,2));
  hg  = hg(:,1); % Gaussian level
 
end	

dt      = t(2)-t(1);
% Calculating covariances
%~~~~~~~~~~~~~~~~~~~~~~~~
R = spec2cov2(S,nr,Ntime-1,dt);

%callFortran = options.method<=0;
callFortran = 0;
if callFortran, % call fortran program
  dens = cov2thpdfexe(R,dt,u,defnr,Nstart,hg,options);
  err  = repmat(nan,size(dens));
else
  [dens,err,terr,options] = cov2thpdf(R,dt,u,defnr,Nstart,hg,options);  
end

def1 = upper(def(1));

if abs(defnr)>1
  f      = createpdf(2);
  f.err  = reshape(err/A,Nx,Ntime);
  f.f    = reshape(dens/A,Nx,Ntime).*der(:,ones(1,Ntime));
  f.x{2} = h;
  f.labx{2}='amplitude [m]';
else
  f   = createpdf(1);
  f.f = dens(:)/A;  
  f.err = err(:)/A;
end
if isfield(spec,'note')
  f.note = spec.note;
end
dt = dt*A;

switch defnr
  case   1, Htxt = sprintf('Density of %sc',def1);
  case  -1, Htxt = sprintf('Density of %st',def1);
  case   2, Htxt = sprintf('Joint density of (%sc,Ac)',def1);
  case  -2, Htxt = sprintf('Joint density of (%st,At)',def1);
  case   3,
    if (t0==tn)
      Htxt = ['Joint density of (',def1,'cf,Ac,',def1,'c=' num2str(tn) ')'];
      f.f   = f.f/dt;
      f.err = f.err/dt;
    else
      Htxt = ['Joint density of (',def1,'cf,Ac)'];
    end
  case  -3,
    if t0==tn
       Htxt = ['Joint density of (',def1,'cb,At,',def1,'t=' num2str(tn) ')'];
       f.f   = f.f/dt;
       f.err = f.err/dt;
    else
      Htxt=['Joint density of (',def1,'cb,At)'];
    end
  case   4, Htxt = ['Joint density of (min(',def1,'cf,',def1,'cb),Ac)'];
  case  -4, Htxt = ['Joint density of (min(',def1,'tf.',def1,'tb),At)'];  
end 

Htxt = [Htxt,sprintf('_{v =%2.5g}',utc)];

f.title = Htxt;
if strcmpi('t',def(1))
  f.labx{1} = 'period [s]';
else
  f.labx{1} = 'wave length [m]';
end

f.x{1}  = t.'*A;
f.options = options;
f.u     = utc;
f.elapsedTime = etime(clock,startTime);

if abs(defnr)>1, 
  try
   [f.cl,f.pl] = qlevels(f.f,[10, 30, 50, 70, 90, 95, 99],f.x{1},f.x{2});
  catch
    warning('WAFO:SPEC2THPDF','Unable to calculate contourlevels')
  end
end

if plotflag  
  pdfplot(f)
end
return % Main
    
function dens = cov2thpdfexe(R,dt,u,defnr,Nstart,hg,options)
  
  Nx    = max(1,length(hg));
  %nr    = size(R,2)-1;
  Ntime = size(R,1);
  
  
  
   rateLHD=3;
   XSPLT = options.xsplit;
   nit   = options.nit;
   speed = options.speed;
   seed  = options.seed;
   SCIS  = abs(options.method); % method<=0
   
   disp('writing data')
   
  filenames0 = writecov(R);
  filenames = {'h.in','reflev.in','dens.out'};
  cleanup(filenames{:})
 
  if ~isempty(hg)
    fid=fopen('h.in','wt');
    fprintf(fid,'%12.10f\n',hg);
    fclose(fid)
  end
   fid=fopen('reflev.in','wt');
   fprintf(fid,'%12.10E \n',u);
   fprintf(fid,'%2.0f \n',defnr);
   fprintf(fid,'%2.0f \n',Ntime);
   fprintf(fid,'%2.0f \n',Nstart);
   fprintf(fid,'%2.0f \n',nit);
   fprintf(fid,'%2.0f \n',speed);
   fprintf(fid,'%2.0f \n',SCIS);
   fprintf(fid,'%2.0f \n',seed);  % select a random seed for rind 
   fprintf(fid,'%2.0f \n',Nx);
   fprintf(fid,'%12.10E \n',dt);
   fprintf(fid,'%2.0f \n',rateLHD);
   fprintf(fid,'%12.10E \n',XSPLT);
   fclose(fid); 
   disp('   Starting Fortran executable.')
   %dos([ wafoexepath 'sp2thpdf70prof.exe'])
   dos([ wafoexepath 'cov2thpdf.exe']);
   %dos([ wafoexepath 'sp2thpdfalan.exe']);   % using EXINV
   %dos([ wafoexepath 'sp2thpdfalan3.exe']);  % using EXINV
   %dos([ wafoexepath 'sp2thpdfalanEX.exe']); % using EXINV
   %dos([ wafoexepath 'sp2thpdfalanFI.exe']); % using FIINV

   %dos([ wafoexepath 'sp2thpdfalanOld.exe']);

   dens = load('dens.out');
   cleanup(filenames0{:},filenames{:})
   return

   function  [pdf,err,terr,options] = cov2thpdf(R,dt,u,def,Nstart,h,options)
%COV2THPDF Joint density of amplitude and period/wave-length
%
% CALL [pdf, err, options] = cov2thpdf(R,dt,u,def,Nstart,h,options)
%
% pdf     = calculated pdf size Nx x Ntime
% err     = error estimate   
% options = requested and actual rindoptions used in integration.
% R       = [R0,R1,R2,R3,R4] column vectors with autocovariance and its
%          derivatives, i.e., Ri (i=1:4) are vectors with the 1'st to 4'th
%          derivatives of R0.  size Ntime x Nr+1
% dt      = time spacing between covariance samples, i.e., 
%           between R0(1),R0(2).
% u       = crossing level
% def     = integer defining pdf calculated:
%           1,  gives half wave period, Tc (default).
%          -1,  gives half wave period, Tt.
%           2,  gives half wave period and wave crest amplitude (Tc,Ac).
%          -2,  gives half wave period and wave trough amplitude (Tt,At).
%           3,  gives crest front period and wave crest amplitude (Tcf,Ac).
%          -3,  gives trough back period and wave trough amplitude (Ttb,At).
%           4,  gives minimum of crest front/back period and wave crest 
%                                           amplitude (min(Tcf,Tcb),Ac).
%          -4,  gives minimum of trough front/back period and wave trough 
%                                           amplitude (min(Ttf,Ttb),At).
% Nstart  = index to where to start calculation, i.e., t0 = t(Nstart)
% h       = vector of amplitudes length Nx or 0
% options = rind options structure defining the integration parameters  
%
% COV2THPDF program computes density of S_i,Hi,T_i in a gaussian process,
%  i.e., quart wavelength (up-crossing to crest) and crest amplitude.
%
% See also  rind, rindoptset  
  
%History:
% Revised pab 22Nov2003  
%  -new inputarguments  
% revised Per A. Brodtkorb 03.03.2003
% -converted from fortran 90 to matlab
% -introduced XcScale instead of CC
% revised Per A. Brodtkorb 04.04.2000
%   - 
% revised Per A. Brodtkorb 23.11.99
%   - fixed a bug in calculating pdf for def = +/- 4
% revised Per A. Brodtkorb 03.11.99
%   - added def = +/-4 
% revised Per A. Brodtkorb 23.09.99
%   - minor changes to covinput
%   - removed the calculation of the transformation to spec2thpdf.m
% by Igor Rychlik
R0 = R(:,1);
R1 = R(:,2);
R2 = R(:,3);
nr = size(R,2)-1;
if nr <3
  R3 = [];
  R4 = [];
else
  R3 = R(:,4);
  R4 = R(:,5);
end

Ntime = length(R0);
Nx    = max(1,length(h));

%C ***** The bound 'infinity' is set to 100*sigma *****
XdInf = 100.d0*sqrt(-R2(1));
XtInf = 100.d0*sqrt(R0(1));
% normalizing constant
%CC=TWPI*SQRT(-R0(1)/R2(1))*exp(u*u/(2.d0*R0(1)) );
XcScale = log(2*pi*sqrt(-R0(1)/R2(1))) + u*u/(2.d0*R0(1));

options.('xcscale') = XcScale;

pdf = zeros(Nx,Ntime);
err = zeros(Nx,Ntime);
terr = err;      
%fxind = zeros(Nx,1);
%err0  = fxind;

%ABSEPS = options.abseps;

opt0 = struct2cell(options);
%opt0 = opt0(1:10);
%opt0 =  {SCIS,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS};
  if (abs(def)<2) 
    NI = 4;
    Nd = 2;
    Nc = 2;
    Mb = 1;
    Nx = 1;
    BIG  = zeros(Ntime+Nc,Ntime+Nc);
    xc   = zeros(Nc,Nx);
    ex   = zeros(1,Ntime+Nc);
    indI = zeros(1,NI);
 %   INFIN = repmat(2,1,NI-1);
    a_up = zeros(Mb,NI-1);
    a_lo = zeros(Mb,NI-1);
    
    
    xc(1,1)=u;
    xc(2,1)=u;
    % INFIN =  INTEGER, array of integration limits flags:  size 1 x Nb   (in)
    %            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
    %            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
    %            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
    %            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
    
    if (def>0) %then
%       INFIN(1:2) = 1;
%       INFIN(3)  = 0;
      a_up(1,1) = u+XtInf ;
      a_lo(1,1)= u;
      a_up(1,2)= XdInf;
      a_lo(1,3)=-XdInf;       
    else
%       INFIN(1:2) = 0;
%       INFIN(3) = 1;
      a_up(1,1)=u ;
      a_lo(1,1)=u-XtInf;
      a_lo(1,2)=-XdInf;
      a_up(1,3)= XdInf ;     
    end %if
    %print *,'Nstart',Nstart
    Nstart=max(2,Nstart)  ;
    %print *,'Nstart',Nstart
    waitTxt =  sprintf('Please wait ...(start at: %s)',datestr(now));   
    h11 = fwaitbar(0,[],waitTxt);

    for Ntd=Nstart:Ntime,
      %CALL COV_INPUT2(BIG,Ntd, R0,R1,R2)
      BIG = covinput(BIG,R0,R1,R2,R3,R4,Ntd,-1); % positive wave period
      Nt  = Ntd-Nd;
      indI(2) = Nt;
      indI(3) = Nt+1;
      indI(4) = Ntd;
      Ntdc    = Ntd+Nc;
      [pdf(Ntd),err(Ntd),terr(Ntd)] = rind(BIG(1:Ntdc,1:Ntdc),ex(1:Ntdc),...
				     a_lo,a_up,indI,xc,Nt,opt0{:});
      waitTxt = sprintf('%s Ready: %d of %d',datestr(now),Ntd,Ntime);
      fwaitbar(Ntd/Ntime,h11,waitTxt);
      
      %disp('sp2thpdf hit a button to continue')
      %pause
    end
    close(h11)
    %err = err + terr;
  else
    XddInf = 100.d0*sqrt(R4(1));
    NI=7; 
    Nd=3;
    Nc=4;
    Mb=2;
    BIG = zeros(Ntime+Nc+1,Ntime+Nc+1);
    xc  = zeros(Nc,Nx);
    ex  = zeros(1,Ntime+Nc+1);
    indI = zeros(1,NI);
 %   INFIN = repmat(2,1,NI-1);
    a_up = zeros(Mb,NI-1);
    a_lo = zeros(Mb,NI-1);
      
    xc(1,1:Nx) = h(1:Nx).';
    xc(2,1:Nx) = u;
    xc(3,1:Nx) = u;
    xc(4,1:Nx) = 0.d0;

    % INFIN =  INTEGER, array of integration limits flags:  size 1 x Nb (in)
    %            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
    %            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
    %            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
    %            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].      
    if (def>0) % then
%       INFIN(2)  = -1;
%       INFIN(4)  = 0;
%       INFIN(5)  = 1;
%       INFIN(6)  = 0;
      a_up(2,1) = 1.d0;   %*h 
      a_lo(1,1) = u;
      a_up(1,2) =  XtInf;   % X(ts) is redundant  
      a_lo(1,2) = -XtInf; 
      a_up(2,2) = 1.d0;   % *h
      a_lo(2,2) = 1.d0;   % *h
      a_up(2,3) = 1.d0;   %*h 
      a_lo(1,3) = u;

      a_lo(1,4) = -XddInf;
      a_up(1,5) =  XdInf;
      a_lo(1,6) = -XdInf ;      
    else %def<0
%       INFIN(2) = -1;
%       INFIN(4) = 1;
%       INFIN(5) = 0;
%       INFIN(6) = 1;
      a_up(1,1) = u ;  
      a_lo(2,1) = 1.d0; %*h
      a_up(1,2) =  XtInf;   % X(ts) is redundant  
      a_lo(1,2) = -XtInf; 
      a_up(2,2) = 1.d0;          % *h
      a_lo(2,2) = 1.d0;          % *h
      a_up(1,3) = u;   
      a_lo(2,3) = 1.d0; %*h
      a_up(1,4) = XddInf;
      a_lo(1,5) = -XdInf;
      a_up(1,6) = XdInf;
    end %if 
    %EPSOLD = ABSEPS;
    Nstart = max(Nstart,3);
    
    waitTxt = sprintf('Please wait ...(start at: %s)',datestr(now));
    h11 = fwaitbar(0,[],waitTxt);

    for tn = Nstart:Ntime,
      Ntd  = tn+1;
      Nt   = Ntd-Nd;
      Ntdc = Ntd+Nc;
      indI(4) = Nt;
      indI(5) = Nt+1;
      indI(6) = Nt+2;
      indI(7) = Ntd;
        
      %opt0{3} =min(sqrt(tn)*EPSOLD*0.5D0,0.1D0);
        
      for ts=2:floor((tn+1)/2.d0),
        %print *,'ts,tn' ,ts,tn
        BIG(1:Ntdc,1:Ntdc) = covinput(BIG(1:Ntdc,1:Ntdc),R0,R1,R2,R3,R4,tn,ts); % positive wave period
        indI(2) = ts-2;
        indI(3) = ts-1;
        [fxind,err0,terr0] = rind(BIG(1:Ntdc,1:Ntdc),ex(1:Ntdc),a_lo, ...
          a_up,indI,xc,Nt,opt0{:});
	  
%         if ts == 4 && tn == 22 || any(isnan(fxind))
%           disp('sp2thpdf hit a button to continue')
%           %pause
%         end
        
	
        switch abs(def)
        case {2}
            % 2,  gives half wave period and wave crest amplitude (Tc,Ac).
          % -2,  gives half wave period and wave trough amplitude (Tt,At).
          if (ts == tn-ts+1) % then
            pdf(1:Nx,tn) = pdf(1:Nx,tn) + fxind.'*dt;
            err(1:Nx,tn) =  err(1:Nx,tn) + (err0(:)*dt).^2;
            terr(1:Nx,tn) =  terr(1:Nx,tn) + terr0(:)*dt;
          else
            pdf(1:Nx,tn) = pdf(1:Nx,tn) + fxind.'*2.d0*dt;
            err(1:Nx,tn) =  err(1:Nx,tn) + (err0(:)*2.0*dt).^2;
            terr(1:Nx,tn) =  terr(1:Nx,tn) + terr0(:)*2.0*dt;
          end %if
        case {3}
          % 3,  gives crest front period and wave crest amplitude (Tcf,Ac).
          % -3,  gives trough back period and wave trough amplitude (Ttb,At).
          pdf(1:Nx,ts) = pdf(1:Nx,ts) + fxind.'*dt;
          err(1:Nx,ts) = err(1:Nx,ts) + (err0(:)*dt).^2;
          terr(1:Nx,ts) =  terr(1:Nx,ts) + terr0(:)*dt;
          if ((ts<tn-ts+1))
            % exploiting the symmetry
            pdf(1:Nx,tn-ts+1) = pdf(1:Nx,tn-ts+1) + fxind.'*dt;
            err(1:Nx,tn-ts+1) = err(1:Nx,tn-ts+1) + (err0(:)*dt).^2;
            terr(1:Nx,tn-ts+1) =  terr(1:Nx,tn-ts+1) + terr0(:)*dt;
          end %if
        case {4}
          % 4,  gives minimum of crest front/back period and wave crest
          %     amplitude (min(Tcf,Tcb),Ac).
          % -4, gives minimum of trough front/back period and wave
          %     trough amplitude (min(Ttf,Ttb),At).
	    
          if (ts == tn-ts+1) % then
            pdf(1:Nx,ts) = pdf(1:Nx,ts)+fxind.'*dt;
            err(1:Nx,ts) = err(1:Nx,ts)+(err0(:)*dt).^2;
            terr(1:Nx,ts) = terr(1:Nx,ts) + terr0(:)*dt;
          else
            pdf(1:Nx,ts) = pdf(1:Nx,ts)+fxind.'*2.0*dt;
            err(1:Nx,ts) = err(1:Nx,ts)+(err0(:)*2.0*dt).^2;
            terr(1:Nx,ts) = terr(1:Nx,ts) + terr0(:)*2.0*dt;
          end %if
        end %select
      end %do                   % ts
      %print *,'Ready: ',tn,' of ',Ntime, ' ABSEPS = ', ABSEPS
      waitTxt = sprintf('%s Ready: %d of %d',datestr(now),tn,Ntime);
      fwaitbar(tn/Ntime,h11,waitTxt);
    end %do                     %tn
    close(h11)
    % err0 is given as 3 standarddeviations of the variability of fxind.
    % thus err is given as a variance
    err = sqrt(err); % convert to standarddeviations
  end %if
  return % main

  

function BIG = covinput(BIG, R0,R1,R2,R3,R4,tn,ts)
%COVINPUT Sets up the covariance matrix  
%
% CALL BIG = covinput(BIG, R0,R1,R2,R3,R4,tn,ts)
%   
%  BIG = covariance matrix for X = [Xt,Xd,Xc] in spec2thpdf problems.
%    
% The order of the variables in the covariance matrix
% are organized as follows: 
% For ts>1:
% ||X(t2)..X(ts),..X(tn-1)||X''(ts) X'(t1) X'(tn)||X(ts) X(t1) X(tn) X'(ts)|| 
% = [Xt                          Xd                    Xc]
%
% For ts<=1: 
% ||X(t2)..,..X(tn-1)||X'(t1) X'(tn)||X(t1) X(tn)||  
% = [Xt                    Xd            Xc]

% where 
%
% Xt = time points in the indicator function
% Xd = derivatives
% Xc = variables to condition on
  
% Computations of all covariances follows simple rules: Cov(X(t),X(s)) = r(t,s),
% then  Cov(X'(t),X(s))=dr(t,s)/dt.  Now for stationary X(t) we have
% a function r(tau) such that Cov(X(t),X(s))=r(s-t) (or r(t-s) will give the same result).
%
% Consequently  Cov(X'(t),X(s))    = -r'(s-t)    = -sign(s-t)*r'(|s-t|)
%               Cov(X'(t),X'(s))   = -r''(s-t)   = -r''(|s-t|)
%               Cov(X''(t),X'(s))  =  r'''(s-t)  = sign(s-t)*r'''(|s-t|)
%               Cov(X''(t),X(s))   =  r''(s-t)   = r''(|s-t|)
% Cov(X''(t),X''(s)) =  r''''(s-t) = r''''(|s-t|)
  
if nargin<8||isempty(ts)
   ts = -1; 
end
if (ts<1) % THEN
   Ntd1 = tn;
 %  N    = Ntd1+2;
   shft = 0;  % def=1 want only crest period Tc
else
   Ntd1 = tn+1;
 %  N    = Ntd1+4;
   shft = 1;  % def=2 or 3 want Tc Ac or Tcf, Ac
end %if
if tn>2
   %do 
   i=1:tn-2;
   %cov(Xt)
   BIG(i,i) = toeplitz(R0(i));% cov(X(ti+1),X(tj+1));

   %cov(Xt,Xc)	
   BIG(i      ,Ntd1+1+shft) = R0(i+1);         %cov(X(ti+1),X(t1))  
   BIG(tn-1-i ,Ntd1+2+shft) = R0(i+1);         %cov(X(t.. ),X(tn))  
   
   %Cov(Xt,Xd)=cov(X(ti+1),x(tj)
   BIG(i,Ntd1-1)    =-R1(i+1);         %cov(X(ti+1),X' (t1))  
   BIG(tn-1-i,Ntd1) = R1(i+1);         %cov(X(ti+1),X' (tn)) 
   %enddo
end
%call echo(big(1:tn,1:tn),tn)
%cov(Xd)
BIG(Ntd1  ,Ntd1  ) = -R2(1);
BIG(Ntd1-1,Ntd1  ) = -R2(tn);     %cov(X'(t1),X'(tn))
BIG(Ntd1-1,Ntd1-1) = -R2(1);

%cov(Xc)
BIG(Ntd1+1+shft,Ntd1+1+shft) = R0(1);        % cov(X(t1),X (t1)) 
BIG(Ntd1+1+shft,Ntd1+2+shft) = R0(tn);      % cov(X(t1),X (tn))
BIG(Ntd1+2+shft,Ntd1+2+shft) = R0(1);        % cov(X(tn),X (tn))
%cov(Xd,Xc)
BIG(Ntd1  ,Ntd1+1+shft) = R1(tn);       %cov(X'(tn),X(t1))     
BIG(Ntd1  ,Ntd1+2+shft) = 0.d0;         %cov(X'(tn),X(tn))
BIG(Ntd1-1,Ntd1+1+shft) = 0.d0 ;        %cov(X'(t1),X(t1))
BIG(Ntd1-1,Ntd1+2+shft) =-R1(tn);      %cov(X'(t1),X(tn))


if (ts>1) % then 

 
  %cov(Xc)
  BIG(Ntd1+1,Ntd1+1) = R0(1);        % cov(X(ts),X (ts)
  BIG(Ntd1+1,Ntd1+2) = R0(ts);       % cov(X(ts),X (t1))
  BIG(Ntd1+1,Ntd1+3) = R0(tn+1-ts);  % cov(X(ts),X (tn))
  BIG(Ntd1+1,Ntd1+4) = 0.d0;         % cov(X(ts),X'(ts))
  
  BIG(Ntd1+2,Ntd1+4) = R1(ts);       % cov(X(t1),X'(ts))
  BIG(Ntd1+3,Ntd1+4) = -R1(tn+1-ts);  %cov(X(tn),X'(ts))
  BIG(Ntd1+4,Ntd1+4) = -R2(1);       % cov(X'(ts),X'(ts))
 
  %cov(Xd)
  BIG(Ntd1-2,Ntd1-1) = -R3(ts);      %cov(X''(ts),X'(t1))
  BIG(Ntd1-2,Ntd1-2) = R4(1);
  BIG(Ntd1-2,Ntd1  ) = R3(tn+1-ts);   %cov(X''(ts),X'(tn))
  %cov(Xd,Xc)
  BIG(Ntd1  ,Ntd1+4) =-R2(tn+1-ts);   %cov(X'(tn),X'(ts))  
  BIG(Ntd1  ,Ntd1+1) = R1(tn+1-ts);   %cov(X'(tn),X (ts))  

  BIG(Ntd1-1,Ntd1+4) =-R2(ts);       %cov(X'(t1),X'(ts))     
  BIG(Ntd1-1,Ntd1+1) =-R1(ts);       %cov(X'(t1),X (ts))  
 
  BIG(Ntd1-2,Ntd1+1) = R2(1);        %cov(X''(ts),X (ts)
  BIG(Ntd1-2,Ntd1+2) = R2(ts);       %cov(X''(ts),X (t1))
  BIG(Ntd1-2,Ntd1+3) = R2(tn+1-ts);   %cov(X''(ts),X (tn))
  BIG(Ntd1-2,Ntd1+4) = 0.d0;         %cov(X''(ts),X'(ts))
  %cov(Xt,Xc)
  if tn>2
    %do 
    i   = 1:tn-2;
    j   = (i+1-ts).';
    tau = abs(j)+1;
    BIG(i,Ntd1+1)   = R0(tau); %cov(X(ti+1),X(ts))  
    BIG(i,Ntd1+4)   = -abs(R1(tau)).*sign(R1(tau).*j); %cov(X(ti+1),X'(ts))  % check this
        
    %Cov(Xt,Xd)=cov(X(ti+1),X(ts))       
    BIG(i,Ntd1-2)    = R2(tau); %cov(X(ti+1),X''(ts))  
				%enddo
  end %tn>2
end %if % ts>1

% make lower triangular part equal to upper 
%for j=1:N-1
%   for i=j+1:N
%      BIG(i,j) = BIG(j,i)
%   end %do
%end %do
lp      = find(tril(ones(size(BIG)),-1)); % indices to lower triangular part
BIGT    = BIG.';
BIG(lp) = BIGT(lp);
% 
return

   