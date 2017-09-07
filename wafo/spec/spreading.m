function D = spreading(th,type,th0,data,w,def)
%SPREADING Directional spreading functions
%
% CALL:  D = spreading(th,type,th0,data,w,def);
%        D = spreading(Nt,type,th0,data,w,def);
%
%      D    = spreading function, struct reminding of spectrum struct   
%      th   = vector of direction angles, (default linspace(-pi,pi,Nt))
%      Nt   = scalar defining the length of th. (default 101)
%      type = type of spreading function, see options below (default 'cos2s')
%      th0  = vector or a scalar defining average direction at every frequency.
%             (length 1 or length == length(w)) (default 0)
%      data = vector of spreading parameters, see options below.
%             (default [15   15   0.52 5    -2.5  0    1    inf])
%      w    = frequency vector, if frequency dependent spreading 
%             (default linspace(0,3,257))
%      def  = 0 No frequency dependence of spreading function                
%             1 Mitsuyasu et al., Hasselman et al (JONSWAP experiment)
%               frequency dependent parametrization (default)
%             
% The different types of spreading functions implemented are:
%  type = 'cos2s'  : cos-2s spreading    N(S)*[cos((th-th0)/2)]^(2*S)  (0 < S) 
%         'box'    : Box-car spreading   N(A)*I( -A < th-th0 < A)      (0 < A < pi)
%         'mises'  : von Mises spreading N(K)*exp(K*cos(th-th0))       (0 < K)
%         'poisson': Poisson spreading   N(X)/(1-2*X*cos(th-th0)+X^2)  (0 < X < 1)
%         'sech2'  : sech-2 spreading    N(B)*0.5*B*sech(B*(th-th0))^2 (0 < B)
%         'wnormal': Wrapped Normal      
%                  [1 + 2*sum exp(-(n*D1)^2/2)*cos(n*(th-th0))]/(2*pi)  (0 < D1)
%         (N(.) = normalization factor)       
%         (the first letter is enough for unique identification)  
%
% Here the S-parameter of the COS-2S spreading function is used as a
% measure of spread. All the parameters of the other distributions are
% related to this S-parameter throug the first Fourier coefficient, R1, of the
% directional distribution as follows: 
%         R1 = S/(S+1) or S = R1/(1-R1).
% where 
%         Box-car spreading  : R1 = sin(A)/A
%         Von Mises spreading: R1 = besseli(1,K)/besseli(0,K), 
%         Poisson spreading  : R1 = X
%         sech-2 spreading   : R1 = pi/(2*B*sinh(pi/(2*B))
%         Wrapped Normal     : R1 = exp(-D^2/2)
%
% Use of data vector:
%   frequency dependent spreading:
%    data =  [spa spb wc ma mb wlim],
%           sp = maximum spread 
%           wc = cut over frequency (usually the peak frequency, wp or fp) 
%           m  = shape parameter defining the frequency dependent
%                spreading parameter, S=S(w), where
%                S(w) = sp *(w/wc)^m, with sp=spa, m=ma for wlim(1) <= w/wc < wlim(2)
%                                          sp=spb, m=mb for wlim(2) <= w/wc < wlim(3)
%                     = 0     otherwise 
%  frequency independent spreading:
%    data =  S,  default value: 15 which corresponds to 
%           'cos2s'  : S=15                      
%           'box'    : A=0.62   
%           'sech2'  : B=0.89      
%           'mises'  : K=8.3      
%           'poisson': X=0.94  
%           'wnormal': D=0.36
%
% The 'cos2s' is the most frequently used spreading in engineering practice.
% Apart from the current meter/pressure cell data in WADIC all
% instruments seem to support the 'cos2s' distribution for heavier sea
% states, (Krogstad and Barstow, 1999). For medium sea states
% a spreading function between 'cos2s' and 'poisson' seem appropriate,
% while 'poisson' seems appropriate for swell.
%   For the 'cos2s' Mitsuyasu et al. parameterized SPa = SPb =
% 11.5*(U10/Cp) where Cp = g/wp is the deep water phase speed at wp and
% U10 the wind speed at reference height 10m. Hasselman et al. (1980)
% parameterized  mb = -2.33-1.45*(U10/Cp-1.17).
% Mitsuyasu et al. (1975) showed that SP for wind waves varies from 
% 5 to 30 being a function of dimensionless wind speed.
% However, Goda and Suzuki (1975) proposed SP = 10 for wind waves, SP = 25
% for swell with short decay distance and SP = 75 for long decay distance.
% Compared to experiments Krogstad et al. (1998) found that ma = 5 +/- eps and
% that -1< mb < -3.5. 
% Values given in the litterature:        [spa  spb  wc   ma   mb      wlim(1:3)  ]
%      (Mitsuyasu: spa == spb)  (cos-2s)  [15   15   0.52 5    -2.5  0    1    inf]
%      (Hasselman: spa ~= spb)  (cos-2s)  [6.97 9.77 0.52 4.06 -2.52 0    1    inf]
%  
%  NOTE: - by specifying NaN's in the data vector default values will be used.
%        - if length(data) is shorter than the parameters needed then the 
%          default values are used
%
% Example:% Set  spa = 10,  wc = 0.43 and spb, ma, mb, wlim to their 
%         % default values, respectively: 
% 
%   data = [10, nan, .43]; 
%   D = spreading(51,'cos2s',0,data);
%        % Frequency dependent direction
%   th0 = linspace(0,pi/2,257)';
%   D = spreading(51,'cos2s',th0,data);
%
% See also  mkdspec, plotspec, spec2spec
 
% References
%  Krogstad, H.E. and Barstow, S.F. (1999)
%  "Directional Distributions in Ocean Wave Spectra"
%  Proceedings of the 9th ISOPE Conference, Vol III, pp. 79-86
%
%  Goda, Y. (1999)
%  "Numerical simulation of ocean waves for statistical analysis"
%  Marine Tech. Soc. Journal, Vol. 33, No. 3, pp 5--14 
%
%  Banner, M.L. (1990)
%  "Equilibrium spectra of wind waves."
%  J. Phys. Ocean, Vol 20, pp 966--984
%
% Donelan M.A., Hamilton J, Hui W.H. (1985)
% "Directional spectra of wind generated waves."
% Phil. Trans. Royal Soc. London, Vol A315, pp 387--407
%
% Hasselmann D, Dunckel M, Ewing JA (1980)
% "Directional spectra observed during JONSWAP."
%  J. Phys. Ocean, Vol.10, pp 1264--1280
%
%  Mitsuyasu, H, et al. (1975)
%  "Observation of the directional spectrum of ocean waves using a
%  coverleaf buoy."
%  J. Physical Oceanography, Vol.5, No.4, pp 750--760

% Some of this might be included in help header:
% cos-2s:
% NB! The generally strong frequency dependence in directional spread
% makes it questionable to run load tests of ships and structures with a
% directional spread independent of frequency (Krogstad and Barstow, 1999).

% Parameterization of B
%    def = 2 Donelan et al freq. parametrization for 'sech2'
%    def = 3 Banner freq. parametrization for 'sech2'
%    (spa ~= spb)  (sech-2)  [2.61 2.28 0.52 1.3  -1.3  0.56 0.95 1.6] 


% Tested on: Matlab 5.3
% History:
% revisd pab 17.06.2001
% - added wrapped normal spreading
% revised pab 6 April 2001
%  - added fcof2par
%  - Fixed the normalization of sech2 spreading
% revised by PAB and IR 1 April 2001: Introducing the azymuth as a
% standard parameter in order to avoid rotations of the directions
% theta. The x-axis is always pointing into the principal direction
% as defined in the spreading function D(omega,theta). The actual
% principal direction is defined by means of field D.phi.
% revised es 06.06.2000, commented away: if ((ma==0) & (mb==0)), ...,
%                    hoping that the check is unnecessary
% revised pab 13.06.2000
%  - fixed a serious bug: made sure -pi<= th-th0 <=pi
% revised pab 16.02.2000
%  -fixed default value for Hasselman parametrization 
% revised pab 02.02.2000
%   - Nt or th may be specified + check on th
%   - added frequency dependence for sech-2
%   - th0 as separate input
%   - updated header info
%   - changed check for nargins
%   - added the possibility of nan's in data vector
% Revised by jr 2000.01.25
% - changed check of nargins
% - frequency dependence only for cos-2s
% - updated information
% By es, jr 1999.11.25

%error(nargchk(0,6,nargin))
narginchk(0,6)
% Default values
%~~~~~~~~~~~~~~~
Nt=101;Nw = 257;
if nargin<1||length(th)<2
  if (nargin>0) && (length(th)==1), 
    Nt=abs(th);
  end
  th = linspace(-pi,pi,Nt).';
elseif abs(th(end)-th(1))<2*pi-eps || abs(th(end)-th(1))>2*pi+eps 
  error('theta must cover all angles -pi -> pi')
else
  Nt = length(th);
end
if Nt<40, 
  warning('WAFO:SPREADING','Number of angles is less than 40. Spreading too sparsely sampled!')
end
if nargin<2||isempty(type), type = 'cos2s'; end
if nargin<3||isempty(th0),  th0  = 0; end
if nargin<5||isempty(w),    w    = linspace(0,3,Nw).'; end
if nargin<6||isempty(def),  def  = 1; end

% Determine the default values

data2=[15   15   0.52 5    -2.5  0    1    inf];
spb=[];wc =[]; ma=[]; mb=[];wlim=[];
nd2 = length(data2);
if (nargin>3) && ~isempty(data), 
  nd  = length(data); 
  ind = find(~isnan(data(1:min(nd,nd2))));
  if any(ind) % replace default values with those from input data
    data2(ind)=data(ind);
  end
  %if ((def==1) & strncmp(type,'cos2s',1) & ((nd<2)|isempty(ind)|~any(ind==2))), 
  %  data2(2)=data2(1); 
  %end
end

if (nd2>0),  spa = data2(1); end
if (nd2>1),  spb = data2(2); end
if def~=0
  if (nd2>2),  wc  = data2(3); end
  if (nd2>3),  ma  = data2(4); end
  if (nd2>4),  mb  = data2(5); end
  if (nd2>5), wlim = data2(6:5+3); end
end
%th = mod(th(:)+pi,2*pi)-pi;
% Make sure theta is from -pi to pi
th      = th(:);
phi0    = th(1)+pi; 
th      = th-phi0;
th(end) = pi;
phi0    = mod(phi0,2*pi);

th0 = th0(:).';
Nt0 = length(th0);
w   = w(:);
Nw  = length(w);

if Nt0==Nw, 
  % frequency dependent spreading and/or
  % frequency dependent direction
  TH = mod(th(:,ones(1,Nw))-th0(ones(Nt,1),:)+pi,2*pi)-pi; % make sure -pi<=TH<pi
  if def==0, def =1;end
elseif Nt0~=1,
  error('The length of th0 must equal to 1 or the length of w')
else
  % If ma==0 and mb==0 then frequency independent spreading is wanted
  %if ((nd2>4) & (ma==0) & (mb==0)), def=0;spb=[];wc =[]; ma=[]; mb=[];wlim=[];end
  TH = mod(th-th0+pi,2*pi)-pi; % make sure -pi<=TH<pi
  if def~=0, % frequency dependent spreading
    TH = TH(:,ones(1,Nw));
  end
end
  
% Createspreading
%~~~~~~~~~~~~~~~~
if def>0
  D  = struct('S',[],'w',w,'theta',th,'type','dir','phi',phi0,'note',[]); % frequency dependent
  if ~isempty(wc),wn = w./wc; end % normalized frequency
else
  D  = struct('S',[],'theta',th,'type','dir','phi',phi0,'note',[]); % frequency independent
end


D.th0  = th0;
D.data = [spa spb wc ma mb wlim]; % save the spreading parameters used
if  def==0||isempty(wc), 
  % no frequency dependent spreading, but possible frequency dependent
  % direction  
  s = spa;
else
  k = find(wlim(3)<wn);
  % Mitsuyasu et. al and Hasselman et. al parametrization   of
  % frequency dependent spreading
  s=[ zeros(length(find(wn<=wlim(1))),1) ; ...
	spa*(wn((wlim(1)<wn) & (wn<=wlim(2)))).^ma ;... 
	spb*(wn((wlim(2)<wn) & (wn<=wlim(3)))).^mb ;... 
	zeros(length(k),1) ].';
  if def==2 && any(k),
     % Donelan et. al. parametrization for B in SECH-2 
    %s(k) =  s(max(find(s~=0)));
    s(k) = s(find(s~=0,1,'last'));
  end
  if def==3 && any(k),
     % Banner parametrization  for B in SECH-2 
    s(k) =  10.^(-0.4+0.8393*exp(-0.567*log(wn(k).^2))); 
  end
end

if any(s<0), error('Spreading: SP value must be larger than 0'),end

if strncmp(type,'cos2s',1),
  if isempty(wc),
    S=s;
  else
    S = s(ones(Nt,1),:);
  end
  D.S = gamma(S+1)/2/sqrt(pi)./gamma(S+1/2).*cos(TH/2).^(2*S);
  D.note='Spreading: cos2s'; 
  return
end

r1 = abs(s./(s+1)); % First Fourier coefficient of the directional spreading function.
%plot(s,r1)
if strncmp(type,'poisson',1),  
  if isempty(wc),
    X = r1;
  else
    X = r1(ones(Nt,1),:);
  end
  if any(X>=1), error('POISSON spreading: X value must be less than 1'),end 
  D.S    = (1-X.^2)./(1-(2*cos(TH)-X).*X)/(2*pi);
  D.note = 'Spreading: Poisson';
  return
end
if strncmp(type,'wnormal',1),
  D1 = repmat(inf,size(r1));
  ind = find(r1); % avoid log of 0
  D1(ind) = sqrt(-2*log(r1(ind)));
  %D1 = sqrt(-2*log(r1));
  if any(D1<0|~isreal(D1)), 
    error('WRAPPED NORMAL spreading: D value must be greater than 0.'),
  end 
  D1 = D1.^2/2;
  ix  = (1:Nt-1).';
  ix2 = ix.^2; 
  if ~isempty(wc), 
    D1 = D1(ones(Nt-1,1),:);
  end
  Nd2 = size(D1,2);
  Fcof = [ ones(1,Nd2)/2;...
	exp(-ix2(:,ones(1,Nd2)).*D1)     ]/pi;
  Pcor = [ ones(1,Nd2 ); exp(sqrt(-1)*ix*TH(1,:))]; % correction term to get
    % the correct integration limits
  Fcof = Fcof(1:Nt,:).*conj(Pcor);
  D.S = real(fft(Fcof));
  D.S(D.S<0)=0;         
  D.note = 'Spreading: Wrapped Normal';
  return
end
% Find distribution parameter from first Fourier coefficient.
par = fcof2par(r1,type);


if strncmp(type,'sech2',1), 
  NB = tanh(pi*par); % Normalization factor.
  NB(NB==0) = 1;     % Avoid division by zero
  if isempty(wc),
    B=par; 
  else
    B = par(ones(Nt,1),:); 
    NB = NB(ones(Nt,1),:);
  end
  D.S = 0.5*B.*sech(B.*TH).^2./NB;
  D.note='Spreading: sech2';
  return
end

if strncmp(type,'mises',1),  
  if isempty(wc),
    K = par; 
  else 
    K = par(ones(Nt,1),:);
  end
  D.S=1./(2*pi*besseli(0,K,1)).*exp(K.*(cos(TH)-1));
  D.note='Spreading: von Mises';
  return
end

if strncmp(type,'box',1),
  if any(par>pi), 
    error('BOX-CAR spreading: The A value must be less than pi'),
  end
  if isempty(wc),
    A=par; 
  else
    A = par(ones(Nt,1),:);
  end
  D.S=1./(2.*A).*(-A<=TH & TH<=A); 
  D.note='Spreading: Box-car';
  return
end

error('TYPE of spreading function unknown')



function x=fcof2par(r1,type)
%FCOF2PAR Fourier coefficients to distribution parameter  

mvrs = version;ix=find(mvrs=='.');
mvnr = str2num(mvrs(1:ix(2)-1));

x = zeros(size(r1));

switch type(1),
  case 's' % sech2
    x0 = [linspace(eps,5,513) linspace(5.0001,100)]';
    x0 = interp1(0.5*pi./(x0.*sinh(.5*pi./x0)),x0,r1(:));
    k1 = find(isnan(x0));
    if any(k1),
      x0(k1)= 0;
      x0(k1)=max(x0);
    end
  
    ix0 = find(r1~=0);
    if mvnr>5.2,
      for ix=ix0,
        x(ix) = abs(fzero('(0.5*pi./(sinh(.5*pi./x)))-x.*(P1)',x0(ix),optimset,r1(ix)));
      end
    else
      for ix=ix0,
        x(ix) = abs(fzero('(0.5*pi./(sinh(.5*pi./x)))-x.*(P1)',x0(ix),[],[],r1(ix)));
      end
    end
  case 'm', % mises
    x0 = [linspace(0,10,513) linspace(10.00001,100)]';
    x0 = interp1(besseli(1,x0,1)./besseli(0,x0,1),x0,r1(:));
    k1 = find(isnan(x0));
    if any(k1),
      x0(k1)= 0;
      x0(k1)=max(x0);
    end
  
    if mvnr>5.2,
      for ix=1:length(r1),
        x(ix) = fzero('besseli(1,x,1)./besseli(0,x,1)-P1',x0(ix),optimset,r1(ix));
      end
    else
      for ix=1:length(r1),
        x(ix) = fzero('besseli(1,x,1)./besseli(0,x,1)-P1',x0(ix),[],[],r1(ix));
      end
    end
  case 'b', % Box-char
    % Initial guess
    x0 = linspace(0,pi+0.1,1025).';
    x0 = interp1(sinc(x0/pi),x0,r1(:));
    x  = x0.';
    
    % Newton-Raphson 
    dx = ones(size(r1));
    iy=0;
    max_count=100;
    ix =find(x);
    while (any(ix) && iy<max_count)
      xi=x(ix);
      dx(ix) = (sin(xi) -xi.*r1(ix)) ./(cos(xi)-r1(ix));
      xi = xi - dx(ix);
      % Make sure that the current guess is larger than zero and less than pi
      %x(ix) = xi + 0.1*(dx(ix) - 9*xi).*(xi<=0) +...
	%  0.38*(dx(ix)-6.2*xi +6.2).*(xi>pi) ;
      % Make sure that the current guess is larger than zero.
      x(ix) = xi + 0.5*(dx(ix) - xi).* (xi<=0);
      
      iy=iy+1;
      ix=find((abs(dx) > sqrt(eps)*abs(x))  &  abs(dx) > sqrt(eps)); 
      %disp(['Iteration ',num2str(iy),...
      %'  Number of points left:  ' num2str(length(ix)) ]),
    end
    if iy == max_count, 
      disp('Warning: fcof2par did not converge.');
      str = 'The last step was:  ';
      outstr = sprintf([str,'%13.8f'],dx(ix));
      fprintf(outstr);
    end
    x(x==0)=eps; % Avoid division by zero  
end
  %semilogy(x,abs(x(:)-x0))
  %plot(x,x(:)-x0)
  
  
  
return


    


