function [fxind, err,terr,exTime,options] = rind(BIG,Ex,Blo,Bup,indI,xc,Nt,varargin)
%RIND Computes multivariate normal expectations
% 
% CALL [E, err,terr,exTime,options] = rind(S,m,Blo,Bup,indI,xc,Nt,options);
%
%        E = expectation/density as explained below size 1 x Nx
%      err = estimated sampling error, with 99% confidence level. size 1 x Nx
%     terr = estimated truncation error.
%   exTime = execution time
%        S = Covariance matrix of X=[Xt;Xd;Xc] size Ntdc x Ntdc (Ntdc=Nt+Nd+Nc)
%        m = the expectation of X=[Xt;Xd;Xc]   size N x 1
%  Blo,Bup = Lower and upper barriers used to compute the integration 
%            limits, Hlo and Hup, respectively. size  Mb x Nb 
%     indI = vector of indices to the different barriers in the  
%            indicator function,  length NI, where   NI = Nb+1 
%            (NB! restriction  indI(1)=0, indI(NI)=Nt+Nd )
%            (default indI = 0:Nt+Nd)
%       xc = values to condition on (default xc = zeros(0,1)), size Nc x Nx
%       Nt = size of Xt             (default Nt = Ntdc - Nc)
%  options = rindoptions structure or named parameters with corresponding
%            values, see rindoptset for details
%
% RIND computes multivariate normal expectations, i.e.,
%    E[Jacobian*Indicator|Condition ]*f_{Xc}(xc(:,ix)) 
% where
%     "Indicator" = I{ Hlo(i) < X(i) < Hup(i), i = 1:N_t+N_d }
%     "Jacobian"  = J(X(Nt+1),...,X(Nt+Nd+Nc)), special case is 
%     "Jacobian"  = |X(Nt+1)*...*X(Nt+Nd)|=|Xd(1)*Xd(2)..Xd(Nd)|
%     "condition" = Xc=xc(:,ix),  ix=1,...,Nx.
%     X = [Xt; Xd; Xc], a stochastic vector of Multivariate Gaussian 
%         variables where Xt,Xd and Xc have the length Nt, Nd and Nc,
%         respectively. 
%     (Recommended limitations Nx,Nt<=100, Nd<=6 and Nc<=10)  
%  
% Multivariate probability is computed if Nd = 0.
%            
% If  Mb<Nc+1 then Blo(Mb+1:Nc+1,:) is assumed to be zero.
% The relation to the integration limits Hlo and Hup are as follows
%
%      Hlo(i) = Blo(1,j)+Blo(2:Mb,j).'*xc(1:Mb-1,ix), 
%      Hup(i) = Bup(1,j)+Bup(2:Mb,j).'*xc(1:Mb-1,ix), 
%
%  where i=indI(j-1)+1:indI(j), j=2:NI, ix=1:Nx
%
% NOTE : RIND is only using upper triangular part of covariance matrix, S
% (except for options.method=0).   
%      
% Examples:% A) Compute Prob{Xi<-1.2} for i=1:5 where
%          %    Xi are zero mean Gaussian variables with covariance
%          %     Cov(Xi,Xj) = 0.3 for i~=j and
%          %     Cov(Xi,Xi) = 1   otherwise
%          % B) Compute expectation E( X1^{+}*X2^{+} ) with random
%          %    correlation coefficient,Cov(X1,X2) = rho2.
%
%  N   = 5; 
%  Blo =-inf; Bup=-1.2; indI=[0 N];              % Barriers
%  A = repmat(Blo,1,N); B = repmat(Bup,1,N); % Integration limits
%  m   = zeros(N,1); rho = 0.3; 
%  Sc  =(ones(N)-eye(N))*rho+eye(N);
%  
%  E0  = rind(Sc,m,Blo,Bup,indI)   % exact prob. 0.001946  A)
%  E1  = rind(triu(Sc),m,A,B)  % same as E0 
%
%  m2   = [0 0]; rho2 = rand(1); 
%  Sc2  = [1 rho2; rho2,1]; 
%  Blo2 = 0; Bup2 = inf; indI2 = [0 2];
%  Nt2  = 0;
%  opt2 = rindoptset('method',1);
%  g2   = inline('(x*(pi/2+asin(x))+sqrt(1-x^2))/(2*pi)'); 
%  
%  E2   = g2(rho2)                          % exact value B)
%  E3   = rind(Sc2,m2,Blo2,Bup2,indI2,[],Nt2)   
%  E4   = rind(Sc2,m2,Blo2,Bup2,indI2,[],Nt2,opt2)
%  E5   = rind(Sc2,m2,Blo2,Bup2,indI2,[],Nt2,'abseps', 1e-4)
%
% See also  prbnormnd, prbnormndpc, rindoptset

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

  
%
%
% Tested on:  Matlab 5.3,DIGITAL UNIX Fortran90 compiler
% History:
% -revised pab July 2007
% -separated error into err and terr.
% -call mexrind2007 instead of mexrindalan24
% -revised pab Nov 2004
%  - removed unused code + Now Nc1c2 is passed to rindalan24  
% -revised pab May 2004
% completely rewritten by pab 18.02.2003
% - initdata is now obsolete and replaced with rindoptset and initoptions
% Note: the order of input arguments changed once again. (hopefully for the
% last time)  
% revised pab 18.02.2003
% - replaced old call to rind.exe with a meximplementation of it.
% revised ir 10.2 2001 File closed after call to Fortran (last fclose)
% revised jr 16.2 2001  (i) speed added in the call of this routine
%                      (ii) The first four lines after the declaration of 
%                           global variables added (call of initdata, speed)
% revised ir 6.2.2001 Adapted to MATLAB 6
% revised jr 1.2.2001 The definition of Nb, Mb changed again. 
% revised ir 15.11.2000 A bug fixed in the definition of Nb Mb.
% revised ir 10.11.2000 Compilation with rind60.f
% revised jr 09.11.2000 The call to init_data replaced with initdata
%% revised pab 23.05.2000 replaced the matlab call with a call to the F90
%         program RINDD.exe
% revised by Per A. Brodtkorb 14.05.1999
%       - No limitations on size of the inputs
%       - enabled recursive calls
%       - fixed some bugs
%       - added some additonal checks
%       - updated to fortran90 
% by  Igor Rychlik 29.10.1998 (PROGRAM RIND11 --- Version 1.0)

%default options
options  = struct('method',3,...
		  'xcscale',0,...
		  'abseps',0,...
		  'releps',1e-3,...
		  'coveps',1e-10,...
		  'maxpts',40000,...
		  'minpts',0,...
		  'seed',floor(rand*1e9),...
		  'nit',1000,...
		  'xcutoff',[],...
		  'xsplit',1.5,...
		  'quadno',[] ,...
		  'speed',[],...
		  'nc1c2',2);
% If just 'defaults' passed in, return the default options in g
if ((nargin==1) && (nargout <= 1) &&  isequal(BIG,'defaults')),
  fxind = options;
  return
end
%error(nargchk(4,inf,nargin));
narginchk(4,inf)
Ntdc = size(BIG,1);

if nargin<6 || isempty(xc),
  xc = zeros(0,1); 
end
Nc = size(xc,1);
if nargin<7 || isempty(Nt),
  Nt = Ntdc-Nc; 
end
[Mb, Nb] = size(Blo);

Nd  = Ntdc-Nt-Nc;
Ntd = Nt+Nd;

if nargin<5 || isempty(indI),
  if Nb~=Ntd
    error('Inconsistent size of Blo and Bup')
  end
  indI = 0:Ntd; 
end
%NI  = length(indI);
%method,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS,seed,NIT,xCutOff
Np = nargin-7;
if (Np>0) % handle various formats for options input
  switch lower(class(varargin{1}))
   case {'char','struct'},
    options = rindoptset(options,varargin{:});
   case {'cell'}
    if Np==1
      options = rindoptset(options,varargin{1}{:});
    else
      error('Invalid options')
    end
   case {'double'}
    % Make sure it is compatible with old calls
    %varargin is a cell vector with double values, i.e.,
    %{method,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS,...
    % seed,NIT,xCutOff,XSPLIT,QUADNO,SPEED};
    ind = findNonEmptyCells(varargin);
    if any(ind)
      opt0 = struct2cell(options);
      opt0(ind) = varargin(ind);
      options = cell2struct(opt0,fieldnames(options));
    end
   otherwise
    error('Invalid options')
  end
end

if isempty(options.xcutoff)
  truncError = 0.1* max(options.abseps);
  Nc1c2 = max(1,options.nc1c2);
  options.xcutoff = max(min(abs(invnorm(truncError/(Nc1c2*2))),8.5),1.2);
  %options.abseps  = max(options.abseps- truncError,0);
  %options.releps  = max(options.releps- truncError,0);
end



% Old call
%if nargin<9 | isempty(SCIS), SCIS = 1; end
%if nargin<10 | isempty(XcScale), XcScale = 0; end
%if nargin<11 | isempty(ABSEPS), ABSEPS = 0; end
%if nargin<12 | isempty(RELEPS), RELEPS = 1e-3; end
%if nargin<13 | isempty(COVEPS), COVEPS = 1e-10; end
%if nargin<14 | isempty(MAXPTS), MAXPTS = 40000; end
%if nargin<15 | isempty(MINPTS), MINPTS = 0; end
%if nargin<16 | isempty(seed),     seed = floor(rand*1e10); end
%if nargin<17 | isempty(NIT),       NIT = 1000; end
%if nargin<18 | isempty(xCutOff), 
%  xCutOff = max(min(abs(invnorm(max(RELEPS,ABSEPS)*10^(-1+0*(Nt>10)))),7),1.2);
%end



%   INFIN  = INTEGER, array of integration limits flags:  size 1 x Nb
%            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
%            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
%            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
%            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
infinity = 37;
dev = sqrt(diag(BIG).');  % std
ind = find(indI(2:end));
INFIN = repmat(2,1,length(indI)-1);
INFIN(ind) = 2 - (Bup(1,ind) > infinity*dev(indI(ind+1)))-...
  2*(Blo(1,ind) <-infinity*dev(indI(ind+1)));

Bup(1,ind) = min(Bup(1,ind) , infinity*dev(indI(ind+1)));
Blo(1,ind) = max(Blo(1,ind) ,-infinity*dev(indI(ind+1)));

if options.method>0
  %opt0 = {SCIS,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS,seed,NIT,xCutOff};
  opt0 = struct2cell(options);
  opt0{1} = mod(opt0{1},10);
  
 
  
  t0 = clock;
 
  try
    if options.method>10
      [fxind,err] = mexrindalan22(BIG,Ex,indI,Blo,Bup,INFIN,xc,Nt,...
				  opt0{1:10});
    else
      % rind with removal of infis and stricter c1c2
      [fxind,err,terr] = mexrind2007(BIG,Ex,indI,Blo,Bup,INFIN,xc,Nt,opt0{[1:10,14]});
    end
  catch	  
    error('Compile mexrind2007 again')
  end 
  exTime = etime(clock,t0);
else
 
  if isempty(options.speed)
    optons.speed = 0;
%    options.speed = 4;
%    options = rindoptset(options);
  end
  
  
  
  t0 = clock;
  try % mexcompiled function
    if options.method==0;
      NIT1 = options.nit;
    else
      NIT1 = -abs(options.method);
    end
    speed   = options.speed;
    seed    = options.seed;
    xcscale = options.xcscale;
    
    
    REPS  = options.releps;
    EPSS    = options.abseps;
    EPS2    = options.coveps;
    xCutOff = options.xcutoff;
    N_int   = options.quadno(1);
    MAXPTS = options.maxpts;
    opt0 = {EPSS,EPS2,xCutOff,REPS,N_int,MAXPTS};
%    dEPSS  = mxGetScalar(prhs(12))
%    dEPS2  = mxGetScalar(prhs(13))
%    dXc    = mxGetScalar(prhs(14))
%    dREPS  = mxGetScalar(prhs(15))
%    dNINT  = NINT(mxGetScalar(prhs(16)))
%    dMAXPTS = NINT(mxGetScalar(prhs(17)))
         
   
    fxind = mexrind71(BIG,Ex,xc,Nt,NIT1,speed,indI,Blo,Bup,seed,xcscale,opt0{:});
%    fxind2 = callRindExe(BIG,Ex,indI,Blo,Bup,xc,Nt,options).';    
%    if any(abs(fxind(:)-fxind2(:))>1e-4)
%      fxind2-fxind
%      disp('diff rind')
%    end
  catch
    error('Compile mexrind71 again')
%    fxind = callRindExe(BIG,Ex,indI,Blo,Bup,xc,Nt,options);
  end
  err = repmat(NaN,size(fxind));
  terr = repmat(NaN,size(fxind));
  exTime = etime(clock,t0);
 
  
end

return % main
  
function ind = findNonEmptyCells(cellArray)
%FINDNONEMPTYCELLS Return index to non-empty cells  
  try % matlab 5.3 or higher
    ind = find(~cellfun('isempty',cellArray)).';
  catch
    % Slow 
    n = numel(cellArray);
    ind = zeros(1,n);
    for ix = 1:n
      ind(ix) = isempty(cellArray{ix});
    end
    ind = find(~ind);
  end
  return % findNonEmptyCells

function [fxind,tid] = callRindExe(BIG,Ex,indI,B_lo,B_up,xc,Nt,options)
%CALLRINDEXE Call rind.exe from wafoexepath (slow)
%
% This is kept just in case  
  if length(Nt)>=2
    Nj = Nt(2);
    Nt = Nt(1);
  else
    Nj=0;
  end
  Nj = min(Nj,max(Nt,0)); % make sure Nj<Nt

  Ntdc = size(BIG,1);
  [Nc, Nx]=size(xc);
  [Mb, Nb]=size(B_lo);
  NI   = length(indI);
  Nd   = Ntdc-Nt-Nc;
  %Ntd  = Nt+Nd;
  
  filenames = {'BIG.in','Ex.in','xc.in','indI.in','B_lo.in','B_up.in','sizeinfo.in','rind.out'};
  
  cleanup(filenames{:})
  
  disp('   Writing data.')
  filename='BIG.in';
  fid=fopen(filename,'wt');
  for ix=1:Ntdc
    fprintf(fid,'%12.10f \n',BIG(ix,:));
  end
  fclose(fid);
  
  filename='Ex.in';
  fid=fopen(filename,'wt');
  fprintf(fid,'%12.10f \n',Ex);
  fclose(fid);
  
  filename='xc.in';
  fid=fopen(filename,'wt');
  fprintf(fid,'%12.10f \n',xc);
  fclose(fid);

  filename='indI.in';
  fid=fopen(filename,'wt');
  fprintf(fid,'%2.0f \n',indI);
  fclose(fid);

  filename='B_lo.in';
  fid=fopen(filename,'wt');
  fprintf(fid,'%12.10f \n',B_lo');
  fclose(fid);

  filename='B_up.in';
  fid=fopen(filename,'wt');
  fprintf(fid,'%12.10f \n',B_up');
  fclose(fid);
  if options.speed>0
    speed   = options.speed;
  else
    speed = 0;
  end
  XSPLT   = options.xsplit;
  RELEPS  = options.releps;
  EPSS    = options.abseps;
  EPS2    = options.coveps;
  xCutOff = options.xcutoff;
  NIT     = options.nit;
  SCIS    = abs(options.method);
  seed    = options.seed;
  N_int   = options.quadno;
  rateLHD = 20;
    
  
  fid=fopen('sizeinfo.in','wt');
  fprintf(fid,'%2.0f \n', speed);
  fprintf(fid,'%2.0f \n', Nt);
  fprintf(fid,'%2.0f \n', Nd);
  fprintf(fid,'%2.0f \n', Nc);
  fprintf(fid,'%2.0f \n', NI);
  fprintf(fid,'%2.0f \n', Mb);
  fprintf(fid,'%2.0f \n', Nx);
  fprintf(fid,'%2.0f \n', NIT);
  fprintf(fid,'%2.0f \n', Nj);
  fprintf(fid,'%2.0f \n', seed); 
  fprintf(fid,'%2.0f \n', SCIS);
  fprintf(fid,'%2.0f \n', rateLHD);
  fprintf(fid,'%12.10f \n', XSPLT);
  % The following variables are used if speed>0
  fprintf(fid,'%12.10f \n', EPSS);    % absolute error
  fprintf(fid,'%12.10f \n', EPS2); % 
  fprintf(fid,'%12.10f \n', xCutOff); % truncation value
  fprintf(fid,'%12.10f \n', RELEPS);  % relativ error
  fprintf(fid,'%2.0f \n', N_int(1));
  fprintf(fid,'%2.0f \n', 1); % minQnr
  fprintf(fid,'%2.0f \n', N_int(1)); % Le2Qnr
  fclose(fid);
  
  disp('   Starting Fortran executable.')
  t0 = clock;
  
  dos([wafoexepath, 'rindd.exe']);
  if nargout>1
    tid=etime(clock,t0);
  end
  fxind=load('rind.out');
  if (Nc>0 && options.xcscale~=0), % scale the result
    CC = exp(options.xcscale);
    fxind = fxind*CC;
  end
  
  cleanup(filenames{:})
  return