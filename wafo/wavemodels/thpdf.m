function f = thpdf(h,Hm0,Tp,dim)
%THPDF Marginal wave height, Hd, pdf for Torsethaugen spectra. 
%
%  CALL: f = thpdf(h,Hm0,Tp,dim)
% 
%  f   = pdf evaluated at h.
%  h   = vectors of evaluation points.
%  Hm0 = significant wave height [m].
%  Tp  = Spectral peak period    [s].
%  dim = 'time'  : Hd distribution in time (default)
%        'space' : Hd distribution in space
%
% THPDF approximates the marginal distribution of Hd, i.e.,
% zero-downcrossing wave height, for a Gaussian process with a Torsethaugen
% spectral density. The empirical parameters of the model is fitted by
% least squares to simulated Hd data for 600 classes of Hm0 and
% Tp. Between 50000 and 150000 zero-downcrossing waves were simulated for
% each class of Hm0 and Tp.
% THPDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% f = thpdf(h,Hm0,Tp);
% dt = 0.4; w = linspace(0,2*pi/dt,256)';
% S = torsethaugen(w,[Hm0 Tp]);
% xs = spec2sdat(S,20000,dt); rate=8; method=1;
% [S,H] = dat2steep(xs,rate,method);
% fk = kdebin(H,{'L2',.5,'inc',128});
% plot(h,f), hold on, pdfplot(fk,'r'), hold off
%
% See also  thvpdf


% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.  
  
% History
% revised pab jan2004  
% By pab 20.12.2000



%error(nargchk(3,4,nargin))
narginchk(3,4)
if nargin<4||isempty(dim),
  dim = 'time';% dim='time'->wtweibpdf, dim='space'->wggampdf
end

if any(Hm0>12| Hm0>(Tp-2)*12/11 )
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end
if any(Tp>20|Tp<3 )
  disp('Warning: Tp is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end
Hrms = Hm0/sqrt(2);

[a b c] = thwparfun(Hm0,Tp,dim);
  f = pdfweibmod(h/Hrms,a,b,c)/Hrms;
  
  return
  % old code kept just in case
% if strncmpi(dim,'t',1)    
%   [a b c] = thwparfun(Hm0,Tp,dim);
%   f = pdfweibmod(h/Hrms,a,b,c)/Hrms;
% else 
%   [a b c] = thgparfun(Hm0,Tp,dim);
%   f = pdfgengam(h/Hrms,a,b,c)/Hrms;
% end

%return

% old code: saved just in case
% Weibull distribution parameters as a function of e2 and h2
% Hrms = Hm0/sqrt(2);
% if nargin<1 |isempty(h), h=linspace(0,4*Hrms).'; end
% if nargin>4,
%   h = h*Hrms;
%   if pdef==1
%     f = hggam(h,Hm0,Tp,eps2,1);
%   else
%     f = htweib(h,Hm0,Tp,eps2,1);
%   end
% elseif pdef==1
%   f = hggam(h,Hm0,Tp,eps2);
% else
%   f = htweib(h,Hm0,Tp,eps2); 
% end

% f=f.f;
% return

% function f = hggam(h,Hm0,Tp,eps2,norm)
%   pardef =2;
%   if pardef == 1,
%     if nargin<4|isempty(eps2),
%       w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
%       %w    = linspace(0,1000/Tp,8*1024+1).'; % adaptiv spacing
%       eps2   = spec2char(torsethaugen(w,[Hm0,Tp]),{'eps2'});
%     end
%     if eps2<0.4 
%       disp('Warning: eps2 is less than 0.4. The pdf returned is questionable')
%     elseif eps2>1.3
%       disp('Warning: eps2 is larger than 1.3. The pdf returned is questionable')
%     end 
%   end

% switch pardef
%   case 1,% pdfgengam distribution parameters as a function of eps2
%     % Best fit by smoothing spline forcing a=1, b=1 and c=2 for eps2=0. 
%     % Then approximate the spline with a rational polynomial
    
%     da1 = [0.3579    0.9601   -1.8152    1.0009];
%     da2 = [1.5625   -0.1537   -1.1390    1.0000];
%     db1 = [-1.8555   3.5375   -3.4056    1.0026]; 
%     db2 = [-2.0855   4.1624   -3.6399    1.0000];
%     dc1 = [6.9572   -6.7932    2.1760    1.9974];
%     dc2 = [3.0979   -3.0781    0.5369    1.0000];
%     A0 = polyval(da1,eps2)./polyval(da2,eps2);
%     B0 = polyval(db1,eps2)./polyval(db2,eps2);
%     C0 = polyval(dc1,eps2)./polyval(dc2,eps2);
%   case 2, % wave height distribution in time
%     global THGPAR
%     if isempty(THGPAR)
%       THGPAR = load('thgpar.mat');
%     end
%     % Generalized Gamma  distribution parameters as a function of Tp, Hm0 
%     A00 = THGPAR.A00s;
%     B00 = THGPAR.B00s;
%     C00 = THGPAR.C00s;

%     Tpp  = THGPAR.Tp;
%     Hm00 = THGPAR.Hm0;
%     [E1, H1] = meshgrid(Tpp,Hm00);
%     method = '*cubic';
%     A0 = interp2(E1,H1,A00,Tp,Hm0,method);
%     B0 = interp2(E1,H1,B00,Tp,Hm0,method);
%     C0 = interp2(E1,H1,C00,Tp,Hm0,method);
%   case 3,% Waveheight distribution in space
%      global THSSPAR
%     if isempty(THSSPAR)
%       THSSPAR = load('thsspar.mat');
%     end
%     % Generalized Gamma  distribution parameters as a function of Tp, Hm0 
%     A00 = THSSPAR.A00s;
%     B00 = THSSPAR.B00s;
%     C00 = THSSPAR.C00s;

%     Tpp  = THSSPAR.Tp;
%     Hm00 = THSSPAR.Hm0;
%     [E1, H1] = meshgrid(Tpp,Hm00);
%     method = '*cubic';
%     A0 = interp2(E1,H1,A00,Tp,Hm0,method);
%     B0 = interp2(E1,H1,B00,Tp,Hm0,method);
%     C0 = interp2(E1,H1,C00,Tp,Hm0,method);
% end
% Hrms = Hm0/sqrt(2);
% f.f    = pdfgengam(h/Hrms,A0,B0,C0)/Hrms;
% return

% function f = htweib(h,Hm0,Tp,eps2,norm)
%  pardef =7;
%   if pardef < 7,
%     if nargin<4|isempty(eps2),
%       w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
%       %w    = linspace(0,1000/Tp,8*1024+1).'; % adaptiv spacing
%       S = torsethaugen(w,[Hm0,Tp]);
%       eps2   = spec2char(S,{'eps2'});
%     end
%     if eps2<0.4 
%       disp('Warning: eps2 is less than 0.4. The pdf returned is questionable')
%     elseif eps2>1.3
%       disp('Warning: eps2 is larger than 1.3. The pdf returned is questionable')
%     end 
%   end
% switch pardef,
%   case 1,% Wtweibpdf distribution parameters as a function of eps2
%     % Best fit by smoothing spline forcing a=1, b=2 and c=0 for eps2=0. 
%     brks = [0 .1 .2 .4 ,.6, .8, 1, 1.1 1.2]';
%     coefa = [0                  0   0.02260415153596   0.99807186986167; ...
%    2.19065400617385                  0   0.02260415153596  1.00033228501527; ...
%    4.34015195156053   0.65719620185215   0.03709199393185    1.00478335417504; ...
%   -1.59533089716870   3.26128737278847   0.80983543882910   1.07321081664798; ...
%   -6.81273221810880   2.30408883448726   1.92291068028425   1.35286675214799; ...
%   -3.69498826658975  -1.78355049637802   2.09369407217829   1.77511058383946; ...
%   13.33514485443956  -4.00054345633187   0.94471547491809   2.09294747228728; ...
%                   0                  0   0.54466112928490   2.16074873007021];
	    
%   coefb = [ 0                  0   0.32503235228616   1.99054481866418; ...
%    3.28321899128157                  0   0.32503235228616    2.02304805389280; ...
%    5.67672309005450   0.98496569738447   0.37964649056830    2.05883450811270; ...
%   -5.29907238080822   4.39099955141717   1.43842344537222   2.21957621884217; ....
%   -5.89663569823287   1.21155612293224   2.55893458024211   2.64050831092684; ...
%   -6.21824739906323  -2.32642529600749   2.43691697455115   3.15358438630669; ...
%   20.19124578481806  -6.05737373544542   0.77134599291113   3.49816479018411; ...
%                   0                  0   0.16560861936659   3.53491689790559];
% coefc =[                0                  0   0.04818579357214       -0.00817761487085; ...
%    2.94432030165157                  0   0.04818579357214    -0.00335903551363; ...
%    4.77660844045250   0.88329609049547   0.09917317190900    0.00440386414523; ...
%   -1.24578770271258   3.74926115476697   1.01096301945323   0.09778320967047; ...
%   -7.70868155645400   3.00178853313943   2.36117295703451   0.43997995813009; ...
%   -3.98346578867600  -1.62342040073298   2.70373824808144   0.97061663841094; ...
%   13.37833291312857  -4.01349987393858   1.57933014446810   1.41455974568850; ...
%                   0                  0   1.17798015707424   1.54573609430905];
%   pa = mkpp(brks,coefa);
%   pb = mkpp(brks,coefb);
%   pc = mkpp(brks,coefc);
%   A0 = ppval(pa,eps2);	    
%   B0 = ppval(pb,eps2);	   
%   C0 = ppval(pc,eps2);	   
%   case 2,% Wtweibpdf distribution parameters as a function of eps2
%     % Best fit by smoothing spline forcing a=1, b=2 and c=0 for eps2=0. 
%     % Then approximate the spline with a rational polynomial.
    
%     da1 = [1.8262  -2.1141   1.0034];
%     da2 = [-0.1371  0.1626   1.3189  -2.0016   1.0000];
%     db1 = [2.4360  -1.6331   1.9932];
%     db2 = [-1.3266  4.3755  -4.6795   2.4763  -1.0450  1.0000];
%     dc1 = [0.3162  -0.1647   0.0803  -0.0104];
%     dc2 = [-1.4325  6.3773 -11.1967  10.1490  -4.7401  1.0000];
%     A0 = polyval(da1,eps2)./polyval(da2,eps2);
%     B0 = polyval(db1,eps2)./polyval(db2,eps2);
%     C0 = polyval(dc1,eps2)./polyval(dc2,eps2);
%   case 3, % Wtweibpdf distribution parameters as a function of eps2 
%     % LS fit to data
%     % best fit to torsethaugen for eps2 = 0.4 to 1.2;
%     % and forcing through A0=1 B0=2 C0=0 for eps2==0) 
%     A0 = sqrt(1+2.0094*eps2.^1.8492);
%     B0 = 1+sqrt(1+3.4117*eps2.^1.4972);
%     C0 = 0.9994*eps2.^1.6728; 
%   case 4,% Wtweibpdf distribution parameters as a function of eps2
%     % best fit to torsethaugen eps2 = 0.4 to 1.2;
%     A0 = 0.7315+eps2.^0.9916;
%     B0 = 2.0628+eps2.^1.1927;
%     C0 = 0.9994*eps2.^1.6728;
%   case 5,% Naess (1985)
%     R  = spec2cov(S);
% %    A0 = sqrt((1-min(R.R)/R.R(1))/2);% Naess (1985)
%     A0 = sqrt((1-min(R.R)/R.R(1))/2)+0.03;% Modified approach broadbanded time
% %    A0 = sqrt((1-min(R.R)/R.R(1))/2)+0.1;% Modified approach broadbanded space                                 
					  
%     B0 = 2;
%     C0 = 0;
%   case 6,% Naess (1985) wave height in space
%     R  = spec2cov(spec2spec(specinterp(S,0.55),'k1d'));
%     %A0 = sqrt((1-min(R.R)/R.R(1))/2); % Naess (1985)
%     A0 = sqrt((1-min(R.R)/R.R(1))/2)+0.03; % modified approach    
%     B0 = 2;
%     C0 = 0; 
%   case 7,% wave height distribution in time
%     global THWPAR
%     if isempty(THWPAR)
%       THWPAR = load('thwpar.mat');
%     end
%     % Truncated Weibull  distribution parameters as a function of Tp, Hm0 
%     A00 = THWPAR.A00s;
%     B00 = THWPAR.B00s;
%     C00 = THWPAR.C00s;

%     Tpp  = THWPAR.Tp;
%     Hm00 = THWPAR.Hm0;
%     [E1, H1] = meshgrid(Tpp,Hm00);
%     method = '*cubic';
%     A0 = interp2(E1,H1,A00,Tp,Hm0,method);
%     B0 = interp2(E1,H1,B00,Tp,Hm0,method);
%     C0 = interp2(E1,H1,C00,Tp,Hm0,method);
 
% end

% Hrms = Hm0/sqrt(2);
% f = wtweibpdf(h/Hrms,A0,B0,C0)/Hrms;

% return
