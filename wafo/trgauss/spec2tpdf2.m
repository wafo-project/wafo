function [f] = spec2tpdf2(spec,utc,def,paramt,varargin)
%SPEC2TPDF2 Density of crest/trough- period or length, version 2.     
%  
%  CALL:  F   = spec2tpdf2(S,u,def,paramt,options);
%
%        F    = density structure.
%        S    = spectral density structure
%        u    = reference level (default the most frequently crossed level).
%       def   = 'Tc',    gives half wave period, Tc (default).
%               'Tt',    gives half wave period, Tt
%               'Lc' and 'Lt' ditto for wave length.
%     paramt  = [t0 tn Nt] where t0, tn and Nt is the first value, last value 
%               and the number of points, respectively, for which
%               the density will be computed. paramt= [5 5 51] implies
%               that the density is computed only for T=5 and
%               using 51 equidistant points in the interval [0,5].
%     options = rind-options structure containing optional parameters
%               controlling the performance of the integration.
%               See rindoptset for details.
%       []    = default values are used.
%
% SPEC2TPDF2 calculates pdf of halfperiods  Tc, Tt, Lc or Lt 
% in a stationary Gaussian transform process X(t), 
% where Y(t) = g(X(t)) (Y zero-mean Gaussian with spectrum given in S). 
% The transformation, g, can be estimated using LC2TR,
% DAT2TR, HERMITETR or OCHITR.  
%
% Example: % The density of Tc is computed by:
%    S = jonswap;
%    paramt = [0 10 51];
%    options = spec2tpdf2('defaults');
%    f = spec2tpdf2(S,[],'Tc',paramt,options);
%    pdfplot(f)
%    hold on, 
%    plot(f.x{:}, f.f+f.err,'r',f.x{:}, f.f-f.err) % estimated error bounds
%    hold off  
%
% See also  rind, rindoptset, spec2cov2, specnorm, dat2tr, dat2gaus, perioddef, wavedef

%History: 
% Tested on : matlab 5.3
% revised pab 24.11.2003
%  updated call to rind
%  -removed nit and speed from input and replaced with rindoptions structure  
% revised pab 
% -removed the printing to screen, replaced with a call to fwaitbar  
% -fixed some bugs: default values for speed, nit and plotflag was not
% properly set  
% revised jr 16.02.2000  nit -> NIT in call of rind and f.nit = NIT; 
% revised ir 10.02.2000 adopted to MATLAB6
% revised pab 23.05.2000 new name spec2tpdf2
% revised by Per A. Brodtkorb 20.06.1999
% by I. Rychlik 01.10.1998 with name wave_th1.m


defaultSpeed   = 9;
defaultoptions = rindoptset('speed',defaultSpeed);
  
if ((nargin==1) && (nargout <= 1) &&  isequal(spec,'defaults')),
  f = defaultoptions;
  return
end 
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
startTime = clock;
if nargin<3||isempty(def)
  def='tc';
end
if strncmpi('l',def,1)
  spec=spec2spec(spec,'k1d');
elseif strncmpi('t',def,1)
  spec=spec2spec(spec,'freq');
else
  error('Unknown def')
end
switch lower(def)
   case  {'tc','lc'},    defnr = 1; % 'tc' or 'lc'
   case  {'tt','lt'},    defnr =-1; % 'tt' or 'lt'
  otherwise ,error('Unknown def')
end			    
[S, xl4, L0, L2]=specnorm(spec);

A = sqrt(L0/L2);

if nargin<5
  options = defaultoptions;
else
  options = rindoptset(defaultoptions,varargin{:});
end

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

if nargin<4||isempty(paramt)
 % z2 = u^2/2;
  z  = -sign(defnr)*u/sqrt(2);
  expectedMaxPeriod = 2*ceil(2*pi*A*exp(z)*(0.5+erf(z)/2)) ;
  paramt = [0, expectedMaxPeriod, 51];
end


% px=gaus2dat([0., u;1, 5],g); 
% px=abs(px(2,2)-px(1,2));
%Nx = 1;
% if 1 % nargin<5||isempty(h)
%   h=px;
%   %h0=0.;
% else
%   h=abs(min(h));
%   h0=h;
%   if h0>0.01*sqrt(L0)
%    % Nx=2;
%     h=[h; px];
%   else
%     h=px;
%    % h0=0.;
%   end
% end

%h=reshape(h,length(h),1);
%hg=tranproc(utc+sign(defnr)*h,g);


t0     = paramt(1);
tn     = paramt(2);
Ntime  = paramt(3);
t      = levels([0, tn/A, Ntime]); % normalized times
Nstart = max(1 + round(t0/tn*(Ntime-1)),2); % index to starting point to
                                     % evaluate

                                     
                                     
                                     
dt = t(2)-t(1);
nr = 2;
R  = spec2cov2(S,nr,Ntime-1,dt);

			    
xc   = [u; u ];
indI = zeros(4,1);
Nd   = 2;
Nc   = 2;
XdInf = 100.d0*sqrt(-R(1,3));
XtInf = 100.d0*sqrt(R(1,1));

B_up  = [u+XtInf, XdInf, 0];
B_lo  = [u,    0, -XdInf];
%INFIN = [1 1 0];
BIG   = zeros(Ntime+2,Ntime+2);
ex    = zeros(Ntime+2,1);
%CC    = 2*pi*sqrt(-R(1,1)/R(1,3))*exp(u^2/(2*R(1,1)));
%  XcScale = log(CC)
XcScale = log(2*pi*sqrt(-R(1,1)/R(1,3)))+(u^2/(2*R(1,1)));
options.xcscale = XcScale;

opt0 = struct2cell(options);

f     = createpdf;
f.f   = zeros(Ntime,1);
f.err = f.f;

h11 = fwaitbar(0,[],sprintf('Please wait ...(start at: %s)',datestr(now)));
for pt = Nstart:Ntime
  Nt      = pt-Nd;
  Ntd     = Nt+Nd;
  indI(2) = Nt;
  indI(3) = Nt+1;
  indI(4) = Ntd;
  Ntdc    = Ntd+Nc;
  % positive wave period  
  BIG = covinput(pt,R,BIG); 
  
  [f.f(pt), f.err(pt)]= rind(BIG(1:Ntdc,1:Ntdc), ex(1:Ntdc),...
			    B_lo,B_up,indI,xc,Nt,opt0{:});

  fwaitbar(pt/Ntime,h11,sprintf('%s Ready: %d of %d',datestr(now),pt,Ntime));
end
close(h11)



switch lower(def)
 case  'tc'
  Htxt = 'Density of Tc';
 case  'tt'
  Htxt = 'Density of Tt';
 case  'lc'
  Htxt = 'Density of Lc';
 case  'lt'
  Htxt = 'Density of Lt';
end 

if strcmpi(def(1),'l')
  xtxt = 'wave length [m]';
else
  xtxt = 'period [s]';
end
Htxt = [Htxt,sprintf('_{v =%2.5g}',utc)];
 
f.title   = Htxt;
f.labx{1} = xtxt;
f.x{1}    = t*A;
f.f       = f.f/A;
f.err     = f.err/A;
f.options = options;
f.u       = utc;
f.elapsedTime = etime(clock,startTime);

if nargout==0 
  pdfplot(f)
end

return


function big = covinput(pt,R,big)
%COVINPUT 
%
% CALL  big = covinput(pt,R,big)
% 
% R = [R0,R1,R2] column vectors with autocovariance and its
%          derivatives, i.e., Ri (i=1:2) are vectors with the 1'st and 2'nd
%          derivatives of R0.  size Ntime x 3
%  
% the order of the variables in the covariance matrix
% are organized as follows: 
% For pt>1:
% ||X(t2)..X(ts),..X(tn-1)|| X'(t1) X'(tn)|| X(t1) X(tn) || 
% = [Xt                          Xd                    Xc]
%
%where 
%
% Xt= time points in the indicator function
% Xd= derivatives
% Xc=variables to condition on

% Computations of all covariances follows simple rules: Cov(X(t),X(s))=r(t,s),
% then  Cov(X'(t),X(s))=dr(t,s)/dt.  Now for stationary X(t) we have
% a function r(tau) such that Cov(X(t),X(s))=r(s-t) (or r(t-s) will give the same result).
%
% Consequently  Cov(X'(t),X(s))    = -r'(s-t)    = -sign(s-t)*r'(|s-t|)
%               Cov(X'(t),X'(s))   = -r''(s-t)   = -r''(|s-t|)
%               Cov(X''(t),X'(s))  =  r'''(s-t)  =  sign(s-t)*r'''(|s-t|)
%               Cov(X''(t),X(s))   =  r''(s-t)   =   r''(|s-t|)
%               Cov(X''(t),X''(s)) =  r''''(s-t) = r''''(|s-t|)


%cov(Xd)
Sdd = -toeplitz(R([1, pt],3));      
%cov(Xc)
Scc = toeplitz(R([1 pt],1));    
%cov(Xc,Xd)
Scd = [0 R(pt,2); -R(pt,2) 0];


if pt>2,
  %cov(Xt)
  Stt = toeplitz(R(1:pt-2,1));  % Cov(X(tn),X(ts))  = r(ts-tn)   = r(|ts-tn|)
  %cov(Xc,Xt) 
  Sct = R(2:pt-1,1).';          % Cov(X(tn),X(ts))  = r(ts-tn)   = r(|ts-tn|)
  Sct = [Sct;Sct(end:-1:1)];
  %Cov(Xd,Xt)
  Sdt = -R(2:pt-1,2).';         % Cov(X'(t1),X(ts)) = -r'(ts-t1) = r(|s-t|)
  Sdt = [Sdt;-Sdt(end:-1:1)];
  N   = pt + 2;
  
  big(1:N,1:N) = [Stt,   Sdt.',   Sct.';...
		  Sdt,   Sdd,     Scd.';...
		  Sct,   Scd,     Scc];
 %error('covinput')
else
  N = 4;
  big(1:N,1:N) =[ Sdd,   Scd.';...
		  Scd,   Scc];
end
% if  pt<0,
%   big
%   pause
%   
% end
return

