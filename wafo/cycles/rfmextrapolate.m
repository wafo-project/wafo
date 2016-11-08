function [Fest,Pout,Fextreme,Fsmooth,Fest0] = rfmextrapolate2(F,Pin,plotflag)

%RFMEXTRAPOLATE Extrapolates a rainflow matrix.
%
% CALL: Fest = rfmextrapolate(F);
%       [Fest,Pout] = rfmextrapolate(F,Pin,plotflag);
%       [Fest,Pout,Fextreme,Fsmooth,Fest0] = rfmextrapolate(F,Pin,plotflag)
%
%   F           = Observed rainflow matrix                        [n,n]
%   Pin         = Input parameters. (optional)           [struct array]
%    .method    = Method for extrapolating the LC.  (see cmat2extralc for options)
%                 'gpd,ml'  : Generalized Pareto distribution. (default)
%                 'exp,mld' : Exponential distribution.  (Linear in lin-log)
%    .u_lev     = Lower and upper levels for extrapolation of levelcrossings. 
%                 (manual choice of levels)          u_lev=[i_min i_max]
%    .LCfrac    = Fraction of LC for choosing thresholds for extraoplation.
%                 (automatic choice of u_lev)       (default 0.05)
%    .h         = Bandwidth for smoothing, h=-1 gives automatic choice.
%    .beta      = Damage exponent.
%    .Lim       = Limits where to use extreme RFM.  (manual choice)
%                 Settings:  Lim.range  Lim.min  Lim.max  (see cmatcombine)
%    .LimRelDam = Threshold for relative damage. 
%                 (automatic choice of Lim)         (default 0.95)
%    .param     = Defines discretization.
%    .PL        = Values of countour lines.
%   plotflag    = 0: Don't plot diagnostic plots, (default)
%                 1: Plot final result of estimated limiting RFM,
%                 2: Plot results of each step in the extrapolation.
%
%   Fest        = Extrapolated RFM.                               [n,n]
%   Pout        = Output parameters.                     [struct array]
%   Fextreme    = Extreme RFM.                                    [n,n]
%   Fsmooth     = Smoothed RFM.  (Kernel smoothing)               [n,n]
%   Fest0       = Extrapolated RFM, before applying limits.       [n,n]
%   
% Extrapolates the level crossing spectrum for high and for low levels. 
% Computes the 'extreme RFM', an approximation of the RFM which 
% is good for large cycles (i.e. for low minima and high maxima). 
% Computes the 'smoothed RFM' by using kernel smoothing.
% The final estimate of the 'extrapolated RFM' (or 'limiting RFM')
% is a combination of Fextreme and Fsmooth.
%
% Example:
%   [G,Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1);
%   xD = mctpsim({G Gh},2000);
%   Frfc = dtp2rfm(xD,64,'CS');
%   Fest = rfmextrapolate(Frfc,[],1);
%   Grfc = mctp2rfm({G Gh});
%   cmatplot({Frfc Fest; Grfc G},4);
%
%   close all;
%
% See also  cmat2extralc, lc2rfmextreme, smoothcmat, cmatcombine

% References:
%
%   Johannesson, P., and Thomas, J-.J. (2001): 
%   Extrapolation of Rainflow Matrices. 
%   Extremes, Vol. 4, pp. 241-262.

% Tested  on Matlab  5.3, 6.5
%
% History:
% Created by PJ (P�r Johannesson) 24-Jul-2000
% Updated by PJ 25-Apr-2003
%   Added input parameter Pin.Lim

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(1,3,ni));

if ni<2, Pin=[]; end
if ni<3, plotflag=[]; end

if isempty(plotflag)
  plotflag=0;
end

n = length(F);

%
% Default Parameters
%

%Pout.paramD = [1 n n];
Pout.param = [1 n n];

Pout.method = 'gpd,ml';
Pout.u_lev = [];
Pout.LCfrac = [];

Pout.h = -1;   % Automatic choice of bandwidth for Kernel smoothing

Pout.Lim = [];
Pout.LimRelDam = [];
Pout.beta = 7;

% Plot-parameters
Pout.PL = [10:20:90 99 99.9 99.99 99.999];

% Copy input parameters Pin to Pout
if ~isempty(Pin)
  Fname = fieldnames(Pin);
  for i = 1:length(Fname)
    Pout = setfield(Pout,Fname{i},getfield(Pin,Fname{i}));
  end
end

%
% Default Parameters for extrapolatin of LC & where to use extreme RFM
%

if isempty(Pout.u_lev) && isempty(Pout.LCfrac)
    Pout.LCfrac = 0.05;
end
if isempty(Pout.Lim) && isempty(Pout.LimRelDam)
    Pout.LimRelDam = 0.95;
end


% 
param = Pout.param;
paramD = [1 n n];
beta = Pout.beta;
PL = Pout.PL;

uD = levels(paramD);
u = levels(param);
du=u(2)-u(1);

%
% Extrapolate LC
%

% Find levels for extrapolation
LCfrac = Pout.LCfrac;
u_lev = Pout.u_lev;

if plotflag >=2, figure, end

if isempty(u_lev)
  slut=0;
  F1=F; i_lowPrev=0; i_highPrev=0;
  while ~slut
    lc = cmat2lc(paramD,F1);
    
    methodFindLev=3;
    if methodFindLev == 1
      [M,Mi] = max(lc(:,2)); Mi = Mi(1);
      dlc = [diff([0; lc(1:Mi-1,2)]); 0; flipud(diff([0; flipud(lc(Mi+1:end,2))]))];
      
      Nlc = sum(dlc~=0); % Number of non-zero in dlc = Number of jumps in lc
      Nextr = floor(0.2*Nlc); % Number of values in the tail of the distribution
      if Nextr<3, Nextr=3; end % At least 3 values in the tail of the distribution
      
      I = find(dlc~=0);
      i_low0 = I(Nextr+1);      % Low level
      i_high0 = I(end-Nextr);   % High level
      
    elseif methodFindLev == 2
      
      [M,Mi] = max(lc(:,2)); Mi = Mi(1);
      dlc1 = diff([0; lc(1:Mi-1,2)]);
      dlc2 = diff([0; flipud(lc(Mi+1:end,2))]);
      
      Nlc1 = sum(dlc1~=0); % Number of non-zero in dlc = Number of jumps in lc
      Nlc2 = sum(dlc2~=0); % Number of non-zero in dlc = Number of jumps in lc
      Nextr1 = max(ceil(0.4*Nlc1),3); % Number of values in the tail of the distribution
      Nextr2 = max(ceil(0.4*Nlc2),3); % Number of values in the tail of the distribution
      % At least 3 values in the tail of the distribution
      
      I1 = find(dlc1~=0);
      I2 = find(dlc2~=0);
      i_low0 = I1(Nextr1+1);      % Low level
      i_high0 = n-I2(Nextr2+1)+1;   % High level
      
      I = [I1; n-flipud(I2)+1];
      
    elseif methodFindLev == 3
      
      [M,Mi] = max(lc(:,2)); M = M(1); Mi=Mi(1);
      
      I = find(lc(:,2)>LCfrac*M);
      i_low0 = I(1);      % Low level
      i_high0 = I(end);   % High level
      
    end
    
    [M,Mi] = max(lc(:,2)); M = M(1); Mi=Mi(1);
    dlc = [diff([0; lc(1:Mi-1,2)]); 0; flipud(diff([0; flipud(lc(Mi+1:end,2))]))];
    Nlc = sum(dlc~=0); % Number of non-zero in dlc = Number of jumps in lc
    Nextr = 3; % Minimum number of values in the tail of the distribution
    I = find(dlc~=0);
    i_low = max(i_low0,I(Nextr+1));      % Low level
    i_high = min(i_high0,I(end-Nextr));   % High level
    
    % Remove cycles with  min>i_high  and  max<i_low
    F1 = F;
    for i = i_high+1:n
      F1(i,:) = 0;
    end
    for j = 1:i_low-1
      F1(:,j) = 0;
    end
    
    slut = (i_high==i_highPrev) & (i_low==i_lowPrev);
    i_highPrev=i_high; i_lowPrev=i_low;
    
  end
else
  i_low = u_lev(1); i_high = u_lev(2);
  % Remove cycles with  min>i_high  and  max<i_low
  F1 = F;
  for i = i_high+1:n
    F1(i,:) = 0;
  end
  for j = 1:i_low-1
    F1(:,j) = 0;
  end
  lc = cmat2lc(paramD,F);
end

if plotflag >=2 % Diagnostic plot
  [M,Mi] = max(lc(:,2)); M = M(1); Mi=Mi(1);
  lc1 =lc; lc1(lc(:,2)==0,2) = 0.1;
  stairs(u(lc1(1:Mi,1)),lc1(1:Mi,2)), hold on
  stairs(u(lc1(Mi+1:end,1)-1),lc1(Mi+1:end,2)), hold on
  %stairs(lc(:,1),lc(:,2)+1), hold on
  %plot(lc(:,1),lc(:,2)+1), hold on
  %lc_min = min(lc(:,2));
  %I = find(dlc~=0);
  v = axis; axis([u([1 n]) 0.9 v(4)]);
  %plot(u(I),0.9*ones(1,length(I)),'*')
  plot(u([i_low i_high]),0.9*ones(1,2),'r.'), 
  plot(u([i_low i_low]),[0.9 v(4)],'r'),
  plot(u([i_high i_high]),[0.9 v(4)],'r'),
  if ~isempty(LCfrac), plot(u([1 n]),LCfrac*[M M],'r'), end
  hold off
  set(gca,'Yscale','log')
  title('Choice of threshold for extrapolation')
  xlabel('level'), ylabel('Number of crossings')
end


% Extrapolate LC

u_lev = [i_low i_high];
if plotflag>=2
  plotflag2 = 1; figure
else
  plotflag2 = 0;
end

method = Pout.method;
[lcEst,Est] = cmat2extralc(paramD,F,u_lev,method,plotflag2);

if plotflag >= 2
  figure
  lcH =lcEst.High; lcH(lcH(:,2)==0,2) = 1e-100;
  lcL =lcEst.Low; lcL(lcL(:,2)==0,2) = 1e-100;
  plot(u(lc(:,1)),lc(:,2),'-'), hold on
  plot(u(lcH(:,1)),lcH(:,2),'r')
  plot(u(lcL(:,1)),lcL(:,2),'r') 
  hold off, grid on  
  set(gca,'YScale','log')
  Mlc = 10^(ceil(log10(max(lc(:,2))))); 
  v=axis; axis([u([1 n]) Mlc/1e8 Mlc])  
  title('Extrapolated level crossings')
  xlabel('level, u')
  ylabel('Level crossing intensity, \mu(u)')
end
Pout.u_lev = u_lev;

%
% Compute extreme RFM
%

lc1 = [0 0; lcEst.lc; n+1 0];

Ga1 = lc2rfmextreme(lc1);
Ga=Ga1(2:n+1,2:n+1);
Ga(isnan(Ga)) = 0;
  
if plotflag >= 2
  figure
  [qlGa PL] = qlevels(F,PL,u,u);

  contour(u,fliplr(u),flipud(F'),qlGa,'b'), hold on
  contour(u,fliplr(u),flipud(Ga'),qlGa,'r'), hold off
  title('Extreme RFM (red), compared with observed RFM (blue)')
  xlabel('min')
  ylabel('Max')
  v = axis; axis([min(u) max(u) min(u) max(u)]);
  axis('square')
end

%
% Smooth the RFM
%

h=Pout.h;

if h~=0
  % Find NOsubzero = Number of subdiagonals equal to zero
  i=0;
  while sum(diag(F,i)) == 0
    i=i+1;
  end
  NOsubzero = i-1;
  
  if h==-1 % Automatic choice of bandwidth
    h_norm = smoothcmat_hnorm(F,NOsubzero);
    % This choice is optimal if the sample is from a normal distribution
    % It usualy oversmooths, therefore one should choose a smaller bandwidth.
    
    h=0.5*h_norm; % Don't oversmooth !!!
  end
  
  [Gs] = smoothcmat(F,1,h,NOsubzero);
  
else
  Gs=F;
end

Pout.h = h;

if plotflag >= 2
  figure
  [qlGs PL] = qlevels(Gs,PL,u,u);

  contour(u,fliplr(u),flipud(Ga'),qlGs,'r'), hold on
  contour(u,fliplr(u),flipud(Gs'),qlGs,'b'), hold off
  title('Smooth RFM (blue), compared with Extreme RFM (red)')
  xlabel('min')
  ylabel('Max')
  v = axis; axis([min(u) max(u) min(u) max(u)]);
  axis('square')
end

% 
% Combine 'Extreme RFM' and  RFMkernel
%

Dnorm = cmat2dam(paramD,F,beta)*100; % Suppose the measurement is 1/100 of the total lifetime
ampF = cmat2amp(paramD,cmat2dmat(paramD,F/Dnorm,beta));
ampGa = cmat2amp(paramD,cmat2dmat(paramD,Ga/Dnorm,beta));

ptyp=2;
Drel = zeros(n,1);
IampF = (ampF(:,2)==0);
IampGa = (ampGa(:,2)==0);
I = ~IampF & ~IampGa;
if ptyp < 3
  Drel(I)=ampGa(I,2)./ampF(I,2);
  Drel(IampF)=1e10;
else
  Drel(I)=ampF(I,2)./ampGa(I,2);
  Drel(IampGa)=1e10;
end

Drel(IampF & IampGa)=0; 
%Drel(isnan(Drel))=0; Drel(isinf(Drel))=1e10;

if plotflag >= 2
  figure
  subplot(2,1,1),plot(ampF(:,1)*du,ampF(:,2)),hold on
  plot(ampGa(:,1)*du,ampGa(:,2),'r'),hold off
  v = axis; axis([0 n/2*du v(3:4)])
  title('Damage per amplitude'), xlabel('amplitude'), ylabel('D0, Dextreme')
  subplot(2,1,2),plot(ampF(:,1)*du,Drel), hold on
  if ptyp<3
    plot([0 n/2*du],0.95*[1 1],'r--'), hold off
    %v = axis; axis([0 n/2*du v(3:4)])
    v = axis; axis([0 n/2*du 0 5])
    title('Relative damage per amplitude'), xlabel('amplitude'), ylabel('Drel = Dextreme / D0')
    if ptyp==2, set(gca,'YDir','reverse'), end
  elseif ptyp==3
    plot([0 n/2*du],1/0.95*[1 1],'r--'), hold off
    v = axis; axis([0 n/2*du 0 5])
    title('Relative damage per amplitude'), xlabel('amplitude'), ylabel('Drel = D0 / Dextreme')
  end
  
end

Lim = Pout.Lim;

if ~isfield(Lim,'range')
    % Where is the extreme RFM valid?
    % Choose limits according to damage.
    % At which amplitude does the extreme RFM give at least 95% of damage of F.
    LimRelDam = Pout.LimRelDam;
    I = find((ampF(:,2)~=0) & (ampGa(:,2)~=0));
    k = min(find(Drel(I)>LimRelDam));
    if ~isempty(k)
        Lim.range = I(k); % Valid fo amplitudes >= Lim.range
    else
        Lim.range=n; % Not valid, don't use Ga !!!
    end
end

[dummy,imax] = max(lcEst.lc(:,2));
maxLC = imax(1);
gap=floor(0.02*n)+1;
if ~isfield(Lim,'min')
    Lim.min = maxLC-gap;
end
if ~isfield(Lim,'max')
    Lim.max = maxLC+gap;
end

[GaC0,LimOut] = cmatcombine(Ga,Gs,Lim);

Pout.Lim = Lim;

% Upper and lower bound of load values,
% due to possible finite endpoints of the LC and hence of the RFM Ga

f_max = sum(Ga);                   % distribution of maxima
jj=min(find([0 fliplr(f_max)]>0)); % Find upper endpoint
jmax = n-jj+2;                     % Largest possible maximum

f_min = sum(Ga');          % distribution of minima
ii=min(find([0 f_min]>0)); % Find lower endpoint
imin = ii-1;               % Smallest possible minimum

J=meshgrid(1:n);  % Columns = maximum
I=J';             % Rows    = minimum
K1 =  (I < imin) | (J > jmax); % Where shall the RFM be zero?
GaC = GaC0;
GaC(K1) = 0;

Pout.imin = imin;
Pout.jmax = jmax;

if plotflag >= 1
  if plotflag>=2, figure, end
  [qlGaC PL] = qlevels(GaC,PL,u,u);

  contour(u,fliplr(u),flipud(F'),qlGaC,'b'), hold on
  contour(u,fliplr(u),flipud(GaC'),qlGaC,'r'), 
  %plot(u([1 n-LimOut.range]),u([LimOut.range n]),'g')
  plot(u([1 n-LimOut.range+1]),u([LimOut.range n]),'g')  % Mod PJ
  plot(u([1 n]),u([LimOut.max LimOut.max]),'g')
  plot(u([LimOut.min LimOut.min]),u([1 n]),'g'),hold off
  title('Limiting RFM (red), compared with observed RFM (blue)')
  xlabel('min')
  ylabel('Max')
  v = axis; axis([min(u) max(u) min(u) max(u)]);
  axis('square')
end


% Compare damage matrices
if plotflag >= 2
  figure
  
  DmatGs = cmat2dmat(paramD,Gs,beta);
  DmatGa = cmat2dmat(paramD,Ga,beta);
  DmatGaC = cmat2dmat(paramD,GaC,beta);
  DmatF = cmat2dmat(paramD,F,beta);
  
  DamF = sum(sum(DmatF));
  DrelGa = sum(sum(DmatGa))/DamF;
  DrelGs = sum(sum(DmatGs))/DamF;
  DrelGaC = sum(sum(DmatGaC))/DamF;

  subplot(2,2,1), cmatplot(u,u,DmatF,13), title(['RFM, Drel = 1'])
  subplot(2,2,2), cmatplot(u,u,DmatGa,13), title(['Extreme RFM, Drel = ' num2str(DrelGa)])
  subplot(2,2,3), cmatplot(u,u,DmatGs,13), title(['Smoothed RFM, Drel = ' num2str(DrelGs)])
  subplot(2,2,4), cmatplot(u,u,DmatGaC,13), title(['Extrapolated RFM, Drel = ' num2str(DrelGaC)])

end

%
% The Result!!!
%

Fest = GaC;

if no>2
  Fest0 = GaC0;
  Fextreme = Ga;
  Fsmooth = Gs;
end

