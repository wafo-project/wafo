function [ch,chtext] = dspec2char(S,varargin)
%DSPEC2CHAR Evaluates directional spectral characteristics 
% 
%  CALL:  [ch,chtext] = dspec2char(S,fact);
% 
%        ch = a cell vector of spectral characteristics
%    chtext = a cell vector of strings describing the elements of ch, see example.
%        S  = Directional spectral density struct with angular frequency
%      fact = vector of factor integers or a string or
%             a cellarray of strings, see below.(default [1])
%
%  DSPEC2CHAR assumes S is a directional spectrum S(w,theta)=
%  S(w)*D(theta,w). If input spectrum is of wave number type, it is
%  transformed into a directional spectrum before the calculations.
%  For many of the parameters the Fourier series expansion of D(theta,w) is used:
%                      M-1
% D(theta(i)) = {1 + 2*sum [an*cos(n*theta(i)) + bn*sin(n*theta(i))]}/(2*pi)
%                      n=1
% where C1 = sqrt(a1^2+b1^2)  and  C2 = sqrt(a2^2+b2^2)
%
%  Order of output is same as order in 'factors'.
%  Input vector 'factors' correspondence:
%
% Factors calculated at every frequency (frequency dependent parameters):
%   1 FMdir  = atan2(b1,a1)                  Mean wave direction.
%   2 FPdir  = atan2(b2,a2)/2                Principal wave direction.
%   3 FSpr   = sqrt(2*(1-C1))                Directional Spread of Mdir
%   4 FSkew  = -C2*sin(2*(Pdir-Mdir))/Spr^3  Circular Skewness of Mdir
%   5 FMSpr  = atan2(sqrt((0.5*b1.^2*(1+a2))-(a1*b1*b2)+...
%               (0.5*a1.^2*(1-a2))),C1^2)    Mean spreading angle
%   6 FLcrst = sqrt((1-C2)/(1+C2))           Long-Crestedness parameter
%   7 FS1    = C1/(1-C1)         Cos^{2S} distribution dispersion parameter, S
%   8 FS2    = [1+3*C2+sqrt(1+(14+C2)*C2)]/(2*(1-C2)) Alternative estimate of S
%   9 FD1    = sqrt(-2*log(C1))  Wrapped Normal distribution parameter, D.
%  10 FD2    = sqrt(-log(C2)/2)  Alternative estimate of D, see spreading.m.
%
% Factors calculated at the peak frequency, fp = 1/Tp: 
%  11 TpMdir =     Mean wave direction at the spectral peak
%  12 TpSpr  =     Directional Spread of TpMdir
%  13 TpSkew =     Skewness of TpMdir
%  14 Wdir   = {theta(i) | [y,i]=max(max(S.S,[],2))}   Main wave direction
%
% Factors calculated from spectrally weighted averages of the
%   Fourier coefficients a and b (frequency independent parameters):
%  15 Wdir2  = {theta(i) | D(theta(i))==max(D(theta))}   Main wave direction
%  16 Mdir   = atan2(b1,a1)                  Mean wave direction.
%  17 Pdir   = atan2(b2,a2)/2                Principal wave direction.
%  18 Spr    = sqrt(2*(1-C1))                Directional Spread of Mdir
%  19 Skew   = -C2*sin(2*(Pdir-Mdir))/Spr^3  Circular Skewness of Mdir
%  20 MSpr   = atan2(sqrt((0.5*b1.^2*(1+a2))-(a1*b1*b2)+...
%               (0.5*a1.^2*(1-a2))),C1^2)    Mean spreading angle
%  21 Lcrst  = sqrt((1-C2)/(1-C2))           Long-Crestedness parameter
%  22 S1     = C1/(1-C1)         Cos^{2s} distribution dispersion parameter, S
%  23 S2     = [1+3*C2+sqrt(1+(14+C2)*C2)]/(2*(1-C2)) Alternative estimate of S
%  24 D1     = sqrt(-2*log(C1))  Wrapped Normal distribution parameter, D.
%  25 D2     = sqrt(-log(C2)/2)  Alternative estimate of D, see spreading.m.
%  26 TMdir  = atan2(b,a)        Mean wave direction Tucker's method.
%
% Note: All angles are given in radians.
%   
%  Examples:
%    S      = demospec('dir');
%    [ch,txt] = dspec2char(S,1:26);        % fact a vector of integers
%    assert(txt, {'FMdir', 'FPdir','FSpr','FSkew','FMSpr','FLcrst','FS1',...
%       'FS2','FD1','FD2','TpMdir','TpSpr','TpSkew','Wdir','Wdir2','Mdir',...
%       'Pdir', 'Spr', 'Skew', 'MSpr', 'Lcrst', 'S1','S2','D1','D2','TMdir'})
%
%    ch0 = cell2struct(ch,txt,2);          % Make a structure     
%    plot(S.w,ch0.FMdir)
%    assert(dspec2char(S,'wdir'){1}, 0, 1e-10)  % fact a string
%    assert(all([dspec2char(S,{'mdir','pdir'}){:}]<1e-10)) % fact a cellarray of strings
%    assert(all([dspec2char(S,'mdir','pdir'){:}]<1e-10))  % strings
%    
%    close all
% 
%  See also  spec2char, spec2bw, spec2mom, spreading  

% References:
% Krogstad, H.E., Wolf, J., Thompson, S.P., and Wyatt, L.R. (1999)
% 'Methods for intercomparison of wave measurements'
% Coastal Enginering, Vol. 37, pp. 235--257
%
% Krogstad, H.E. (1982)
% 'On the covariance of the periodogram'
% Journal of time series analysis, Vol. 3, No. 3, pp. 195--207
%
% Tucker, M.J. (1993)
% 'Recommended standard for wave data sampling and near-real-time processing'
% Ocean Engineering, Vol.20, No.5, pp. 459--474
%
% Young, I.R. (1999)
% "Wind generated ocean waves"
% Elsevier Ocean Engineering Book Series, Vol. 2, pp 239
%
% Benoit, M. and Goasguen, G. (1999)
% "Comparative evalutation of directional wave analysis techniques applied
% to field measurements", In Proc. 9'th ISOPE conference, Vol III, pp 87-94.

% Tested on: Matlab 5.2

% History: 
%revised pab 29.06.2001
% - factors 15:25 are now calculated from the scaled Fourier coefficients of
% D(theta) = int S(w,theta) dw instead of  D(theta) = int D(w,theta) dw 
%revised pab 22.06.2001
% - moved code from wspecplottest into dspec2char
% - Fixed bugs: The parameters are now calculated from the correct
%     Scaled Fourier coefficients. S.phi is now taken into account.
% - Added frequency dependent parameters, skewness and kurtosis
%   parameters and dispersion parameter for the wrapped normal spreading function.
%by Vengatesan Venugopal [V.Venugopal@hw.ac.uk] 19.06.2001

% Options not implemented
%  FKurt  = 2*(C2*cos(2*(Pdir-Mdir))-4*C1+3)/Spr^4 Circular kurtosis of FMdir
%  TpKurt =     Kurtosis of TpMdir
%  Kurt   = 2*(C2*cos(2*(Pdir-Mdir))-4*C1+3)/Spr^4 Circular kurtosis of Mdir
%  TSpr   = sqrt(2*(1-C1))    Directional Spread of TMdir  (Wrong?)

transform2degrees = 0;

switch nargin
  case 0, return
  case 1, 
  case 2,
    fact = varargin{1};
    if ischar(fact), fact = {fact}; end
  otherwise
    fact = varargin;
end

tfact ={'FMdir','FPdir','FSpr','FSkew','FMSpr','FLcrst','FS1','FS2','FD1','FD2',...
    'TpMdir','TpSpr','TpSkew','Wdir','Wdir2', ...
'Mdir','Pdir','Spr','Skew','MSpr','Lcrst','S1','S2','D1','D2','TMdir','TSpr'};  



if iscell(fact)
  N     = length(fact(:));
  nfact = zeros(1,N);
  ltfact = char(lower(tfact));
  for ix=1:N,
    ind = strmatch(lower(char(fact(ix))),ltfact,'exact');
    if length(ind)==1,
      nfact(ix)=ind;
    else
      error(['Not a valid factor: ' fact{ix}]);
    end
  end 
else
  nfact = fact;
end
if any(nfact>26 | nfact<1)
  error('Factor outside range (1,...,26)');
end

vari = 'w';
if isfield(S,'k2d')
  S = spec2spec(S,'dir');
elseif isfield(S,'theta') 
  S = ttspec(S,'w','r'); % Make sure it is directional spectrum given in radians
else
  error('Directional spectra required!');
end

phi = 0;
if isfield(S,'phi') && ~isempty(S.phi)
  phi = S.phi;
end



w       = S.(vari);
w       = w(:);
theta   = S.theta(:)-phi;
Dtf     = S.S;
[Nt Nf] = size(Dtf);

if length(w)~=Nf, 
  error('Length of frequency vector S.f or S.w must equal size(S.S,2)');
end
if length(theta)~=Nt, 
  error('Length of angular vector S.theta must equal size(S.S,1)');
end

Sf      = simpson(S.theta,Dtf,1);
ind     = find(Sf);

%Directional distribution  D(theta,w) = S(theta,w)/S(w)
Dtf(:,ind) = Dtf(:,ind)./Sf(ones(Nt,1),ind); 

if 0,
  Dtheta   = simpson(w,Dtf,2); %Directional spreading, D(theta) = int D(w,theta) dw
else
  Dtheta   = simpson(w,S.S,2); %Directional spreading, D(theta) = int S(w,theta) dw
end
Dtheta     = Dtheta/simpson(S.theta,Dtheta);

%[y,ind] = max(Dtheta);
%Wdir    = theta(ind); % main wave direction



% Find Fourier Coefficients of Dtheta and Dtf
M = 3; % No of harmonics-1
[a,b]   = fourier(theta,Dtheta,2*pi,10);
[aa,bb] = fourier(theta,Dtf,2*pi,10);
if 0, % Alternatively
  Fcof = 2*ifft(Dtheta);
  Pcor = [1; exp(sqrt(-1)*(1:M-1).'*theta(1))]; % correction term to get
  % the correct integration limits
  Fcof = Fcof(1:M,:).*Pcor;
  a = real(Fcof(1:M));
  b = imag(Fcof(1:M));
  Fcof = 2*ifft(Dtf);
  Fcof = Fcof(1:M,:).*Pcor(:,ones(1,Nf));
  aa = real(Fcof(1:M,:));
  bb = imag(Fcof(1:M,:));
end

if 0,
  % Checking the components
  P1toM = zeros(Nt,M);
  P1toM(:,1) =  a(1)/2; % mean level
  
  for ix=2:M,
    % M-1 1st harmonic
    P1toM(:,ix) = a(ix)*cos((ix-1)*theta)+b(ix)*sin((ix-1)*theta);
  end
  Psum = sum(P1toM,2);
  subplot(2,1,1);
  plot(theta,Dtheta-Psum,'o','MarkerSize',2), legend('Dteta-Psum');
  subplot(2,1,2);
  plot(theta,Dtheta,'o','MarkerSize',2); hold on
  plot(theta,Psum,'m'); hold off;legend('Dteta','Psum');
  pause 
end   

% The parameters below are calculated for 
a  = pi*a;  b  = pi*b;
aa = pi*aa; bb = pi*bb;


%Fourier coefficients for D(theta)  
%a0=a(1);
a1=a(2);a2=a(3);
b1=b(2); b2=b(3);

%Fourier coefficients for D(theta,w)
%aa0 = aa(1,:); 
aa1 = a(2,:);aa2 = aa(3,:);
bb1 = bb(2,:); bb2 = bb(3,:);

FC1 = sqrt(aa1.^2+bb1.^2);
FC2 = sqrt(aa2.^2+bb2.^2);

%plot(w,FC1,w,FC2),legend('FC1','FC2')


FMdir = atan2(bb1,aa1);                          % Mean wave direction
FPdir = 0.5*atan2(bb2,aa2);                      % Principal wave direction
FSpr  = sqrt(2*abs(1-FC1));                      % Directional spread
FSkew = -FC2.*sin(2*(FPdir-FMdir))./(FSpr.^3);   % Skewness of Mdir
%FKurt = 2*(FC2.*cos(2*(FPdir-FMdir))-4.*FC1+3)./(FSpr.^4); % Kurtosis of Mdir

 
% Mean spreading angle 
FMSpr = atan2(sqrt(0.5*bb1.^2.*(1+aa2)-aa1.*bb1.*bb2+0.5.*aa1.^2.*(1-aa2)),FC1.^2);

% Long-Crestedness parameter
FLcrst  = sqrt((1-FC2)./(1+FC2));

%Estimates of Cos^{2s} distribution dispersion parameter, S 
FS1 = FC1./(1-FC1);
FS2 = (1+3*FC2+sqrt(1+(14+FC2).*FC2))./(2*(1-FC2));
%Estimates of Wrapped Normal distribution parameter, D.
FD1     = repmat(-inf,size(FC1));
ind     = find(FC1); % avoid log(0)
FD1(ind)= sqrt(-2*log(FC1(ind)));  

FD2     = repmat(-inf,size(FC2));
ind     = find(FC2); % avoid log(0)
FD2(ind)= sqrt(-log(FC2(ind))/2); 

[y,ind] = max(max(S.S,[],1)); % Index to spectral peak.

TpMdir = FMdir(ind); % Mean wave direction at the spectral peak
TpSpr  = FSpr(ind);  % Directional Spread of TpMdir
TpSkew = FSkew(ind); % Skewness of TpMdir
%TpKurt = FKurt(ind); % Kurtosis of TpMdir


[y,ind] = max(max(S.S,[],2)); % Index to spectral peak.
Wdir    = mod(theta(ind)+pi,2*pi)-pi; % main wave direction
[y,ind] = max(Dtheta);
Wdir2   = mod(theta(ind)+pi,2*pi)-pi; % main wave direction



C1 = sqrt(a1.^2+b1.^2);
C2 = sqrt(a2.^2+b2.^2);

     

Mdir = atan2(b1,a1);                               % Mean wave direction
Pdir = 0.5*atan2(b2,a2);                           % Principal wave direction
Spr  = sqrt(2*abs(1-C1));                          % Directional spread
Skew = -C2.*sin(2*(Pdir-Mdir))./(Spr.^3);          % Skewness of Mdir
%Kurt = 2*(C2.*cos(2*(Pdir-Mdir))-4.*C1+3)./(Spr.^4); % Kurtosis of Mdir
  
% Mean spreading angle 
MSpr = atan2(sqrt(0.5*b1.^2.*(1+a2)-a1.*b1.*b2+0.5.*a1.^2.*(1-a2)),C1.^2);

% Long-Crestedness parameter
Lcrst  = sqrt((1-C2)./(1+C2));

%Estimates of Cos^{2s} distribution dispersion parameter, S 
S1 = C1./(1-C1);
S2 = (1+3*C2+sqrt(1+(14+C2).*C2))./(2*(1-C2));
%Estimates of Wrapped Normal distribution parameter, D.
D1     = sqrt(abs(2*log(C1)));  
D2     = sqrt(abs(log(C2)/2)); 


% ---------------------------------------------------------------------------   
% TUCKER PROCEDURE for MDIR and UI (pp202)
         

[Sff,ind] = max(S.S(:,:));
%Wwdir     = theta(ind); % main wave direction

%thetam = (ind-1).'*2*pi/(Nt-1)-pi;
thetam = theta(ind);
Sff = Sff';
bot = sum(Sff);
a   = sum(Sff.*cos(thetam))/bot;
b   = sum(Sff.*sin(thetam))/bot;

TMdir = atan2(b,a);
%UI   = sqrt(a.^2+b.^2);
%TSpr = sqrt(2*(1-UI));  %This is not correct  

ch ={FMdir,FPdir,FSpr,FSkew,FMSpr,FLcrst,FS1,FS2,FD1,FD2,...
    TpMdir,TpSpr,TpSkew,Wdir,Wdir2,Mdir,Pdir,Spr,Skew,....
    MSpr,Lcrst,S1,S2,D1,D2,TMdir}; 

% Make sure the angles are between -pi<= theta < pi
ind = [1:3 5 11:12 14:18 20 26];
for ix = ind,
  ch{ix} = mod(ch{ix}+pi,2*pi)-pi;
end 
if transform2degrees, % Change from radians to degrees  
  for ix = ind,
    ch{ix} = ch{ix}*180/pi;
  end 
end

% select the appropriate values
ch     = ch(nfact);
chtext = tfact(nfact);

return

