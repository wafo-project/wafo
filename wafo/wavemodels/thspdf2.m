function [f,varargout] = thspdf2(Hd,Scf,Hm0,Tp,normalizedInput)
%THSPDF2 Joint (Scf,Hd) PDF for linear waves with Torsethaugen spectra. 
%
%  CALL: f = thspdf2(Hd,Scf,Hm0,Tp)
% 
%   f   = pdf struct evaluated at meshgrid(Scf,Hd)
%   Hd  = zero down crossing wave height (vector)
%   Scf = crest front steepness (vector) 
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
%
% THSPDF2 approximates the joint distribution of (Scf, Hd), i.e., crest
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for a Gaussian process with a
% Torsethaugen spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Scf,Hd) data for 600 classes of
% Hm0 and Tp. Between 40000 and 200000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THSPDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2)); 
% s = linspace(0,6*1.25*Hm0/Tp^2);
% f = thspdf2(h,s,Hm0,Tp);
% w = linspace(0,40,5*1024+1).';
% S = torsethaugen(w,[Hm0 Tp]);
% dt = 0.3;
% x = spec2sdat(S,80000,.2); rate = 8;
% [si,hi] = dat2steep(x,rate,2);
% fk = kdebin([si,hi],{'L2',.5,'inc',128});
%  fk.title = f.title; fk.labx = f.labx; 
% plot(si,hi,'.'), hold on
% pdfplot(f),pdfplot(fk,'r'),hold off
%
% See also  thsspdf

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.  

% Adapted to  cssmooth  by GL Feb 2011  
% History
% revised pab 09.08.2003
% changed input  
% validated 20.11.2002
% By pab 20.12.2000

%error(nargchk(4,5,nargin))
narginchk(4,5)
if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end


[V,H] = meshgrid(Scf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = thspdf(H,V,Hm0,Tp,normalizedInput);

 f.x = {Scf(:),Hd(:)};
 
if (normalizedInput)
  f.labx={'Scf', 'Hd'};
  f.norm = 1;
else
  f.norm=0;
  f.labx={'Scf', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Scf) in time';
f.note = ['Torsethaugen Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp)];
[f.cl,f.pl] = qlevels(f.f);

return 
% old call
% if Hm0>11| Hm0>(Tp-2)*12/11 
%   disp('Warning: Hm0 is outside the valid range')
%   disp('The validity of the Joint (Hd,Scf) distribution is questionable')
% end
% if Tp>20|Tp<3 
%   disp('Warning: Tp is outside the valid range')
%   disp('The validity of the Joint (Hd,Scf) distribution is questionable')
% end
% 
% global THSPAR
% if isempty(THSPAR)
%   THSPAR = load('thspar.mat');
% end
% % Gamma distribution parameters as a function of Tp Hm0 and h2
% A11 = THSPAR.A11s;
% B11 = THSPAR.B11s;
% 
% % Waveheight distribution in time
% 
% if 0,
%   % Truncated Weibull  distribution parameters as a function of Tp, Hm0 
%   global THWPAR
%   if isempty(THWPAR)
%     THWPAR = load('thwpar.mat');
%   end
%   A00 = THWPAR.A00s;
%   B00 = THWPAR.B00s;
%   C00 = THWPAR.C00s;
% else
% 
%   % Truncated Weibull  distribution parameters as a function of Tp, Hm0 
%   A00 = THSPAR.A00s;
%   B00 = THSPAR.B00s;
%   C00 = THSPAR.C00s;
% end
% 
% Tpp  = THSPAR.Tp;
% Hm00 = THSPAR.Hm0;
% h2   = THSPAR.h2;
% 
% 
% 
% 
% w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
% %w    = linspace(0,10,2*1024+1).'; 
% S = torsethaugen(w,[Hm0,Tp]);
% ch   = spec2char(S,{'Tm02','eps2'});
% Tm02 = ch(1);
% eps2 = ch(2);
% 
% Hrms = Hm0/sqrt(2);
% Vrms = 1.25*Hm0/(Tm02^2); % Erms
% 
% if nargin<1 |isempty(v), v=linspace(0,4*Vrms); end
% if nargin<2 |isempty(h), h=linspace(0,4*Hrms); end
% 
% if nargin>4,
%   h = h*Hrms;
%   v = v*Vrms;
% end
% 
% %Fh = thpdf(h(:)/Hrms,Hm0,Tp,eps2,1);
% 
% [A0, B0, C0] = thwparfun(Hm0,Tp,'time');
% 
% method = '*cubic';% Faster interpolation
% [E1, H1, H2] = meshgrid(Tpp,Hm00,h2);
% A1 = exp(cssmooth(h2,interp3(E1,H1,H2,log(A11),Tp(ones(size(h2))),...
%     Hm0(ones(size(h2))) ,h2,method),1,h/Hrms,1));
% B1 = exp(cssmooth(h2,interp3(E1,H1,H2,log(B11),Tp(ones(size(h2))),...
%      Hm0(ones(size(h2))),h2,method),1,h/Hrms,1));
% 
% [V1 H1] = meshgrid(v/Vrms,h/Hrms);
% [V1 A1] = meshgrid(v/Vrms,A1);
% [V1 B1] = meshgrid(v/Vrms,B1);
% 
% f = createpdf(2);
% f.title = 'Joint distribution of (Hd,Scf) in time';
% 
% if nargin<5 
%   
%    f.f = wtweibpdf(H1,A0,B0,C0).*pdfgam(V1,A1,B1)/Vrms/Hrms;
%  % f.f = repmat(Fh/Hrms,[1 length(v)]).*pdfgam(V1,A1,B1)/Vrms;
%   f.x = {v,h};
%   f.norm=0;
%   f.labx={'Scf', 'Hd [m]'};
% else
%   f.f = wtweibpdf(H1,A0,B0,C0).*pdfgam(V1,A1,B1);
%   %f.f = repmat(Fh,[1 length(v)]).*pdfgam(V1,A1,B1);
%   f.x = {v/Vrms,h/Hrms};
%   f.labx={'Scf', 'Hd'};
%   f.norm = 1;
% end
% f.f(find(isnan(f.f)|isinf(f.f) ))=0;
% 
% f.note = ['Torsethaugen Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp)];
% [f.cl,f.pl] = qlevels(f.f);
% if nargout>1,
%   fA      = createpdf(2);
%   fA.x    = {Tpp,Hm00};
%   fA.labx = {'Tp', 'Hm0'};
%   fA(3)   = fA(1);
%   fA(2)   = fA(1);
%   
%   fA(1).f    = A00;
%   fA(2).f    = B00;
%   fA(3).f    = C00;
%   
%   fA(1).title = 'wtweibpdf parameter A';
%   fA(2).title = 'wtweibpdf parameter B';
%   fA(3).title = 'wtweibpdf parameter C';
%   
%   txt1 = 'The Wtweibpdf  distribution Parameter ';
%   txt2=[' for wave heigth in time as a function of Tp and Hm0 for' ...
% 	'the Torsethaugen spectrum'];
%   fA(1).note =[txt1 'A' txt2];
%   fA(2).note =[txt1 'B' txt2];
%   fA(3).note =[txt1 'C' txt2];
%   
%   tmp = [A00(:) B00(:) C00(:)];
%   ra  = range(tmp);
%   st  = round(min(tmp)*100)/100;
%   en  = max(tmp);
%   for ix = 1:3,
%     fA(ix).cl   = st(ix):ra(ix)/20:en(ix);
%   end
% end
% if nargout>2,
%   fB      = createpdf(3);
%   fB.x    = {Tpp,Hm00,h2};
%   fB.labx = {'Tp','Hm0', 'h'};
%   fB(2)   = fB(1);
%   
%   fB(1).f = A11;
%   fB(2).f = B11;
%   
%   txt11 = 'The conditonal pdfgam distribution Parameter ';
%   txt22 = [' for Scf given h=Hd/Hrms in time as function of Tp' ...
% 	' and Hm0 for the Torsethaugen spectrum'];
%   fB(1).title = 'pdfgam parameter A';
%   fB(2).title = 'pdfgam parameter B';
%   fB(1).note = [txt11,'A',txt22];
%   fB(2).note = [txt11,'B',txt22];
%   
%   %fB(2).note = ['The conditonal pdfgengam distribution Parameter B(h)/Hrms', ...%	' for crest front steepness as a function of Tp,Hm0 and',...
%   %	' h=Hd/Hrms for the Torsethaugen spectrum'];
%   tmp= [A11(:) B11(:)];
%   ra = range(tmp);
%   st = round(min(tmp)*100)/100;
%   en = max(tmp);
%   for ix=1:2
%     fB(ix).cl   = st(ix):ra(ix)/20:en(ix);
%   end
% end
% return
% 
% 
% 
% 
% 
% 
% % %Old calls
% % %method ='spline';
% % switch 1
% %  case 3,% Best fit by smoothing spline
% %     brks = [0 .1 .2 .4 ,.6, .8, 1, 1.1 1.2]';
% %     coefa = [0                  0   0.02260415153596   0.99807186986167; ...
% %    2.19065400617385                  0   0.02260415153596  1.00033228501527; ...
% %    4.34015195156053   0.65719620185215   0.03709199393185    1.00478335417504; ...
% %   -1.59533089716870   3.26128737278847   0.80983543882910   1.07321081664798; ...
% %   -6.81273221810880   2.30408883448726   1.92291068028425   1.35286675214799; ...
% %   -3.69498826658975  -1.78355049637802   2.09369407217829   1.77511058383946; ...
% %   13.33514485443956  -4.00054345633187   0.94471547491809   2.09294747228728; ...
% %                   0                  0   0.54466112928490   2.16074873007021];
% % 	    
% %   coefb = [ 0                  0   0.32503235228616   1.99054481866418; ...
% %    3.28321899128157                  0   0.32503235228616    2.02304805389280; ...
% %    5.67672309005450   0.98496569738447   0.37964649056830    2.05883450811270; ...
% %   -5.29907238080822   4.39099955141717   1.43842344537222   2.21957621884217; ....
% %   -5.89663569823287   1.21155612293224   2.55893458024211   2.64050831092684; ...
% %   -6.21824739906323  -2.32642529600749   2.43691697455115   3.15358438630669; ...
% %   20.19124578481806  -6.05737373544542   0.77134599291113   3.49816479018411; ...
% %                   0                  0   0.16560861936659   3.53491689790559];
% % coefc =[                0                  0   0.04818579357214       -0.00817761487085; ...
% %    2.94432030165157                  0   0.04818579357214    -0.00335903551363; ...
% %    4.77660844045250   0.88329609049547   0.09917317190900    0.00440386414523; ...
% %   -1.24578770271258   3.74926115476697   1.01096301945323   0.09778320967047; ...
% %   -7.70868155645400   3.00178853313943   2.36117295703451   0.43997995813009; ...
% %   -3.98346578867600  -1.62342040073298   2.70373824808144   0.97061663841094; ...
% %   13.37833291312857  -4.01349987393858   1.57933014446810   1.41455974568850; ...
% %                   0                  0   1.17798015707424   1.54573609430905];
% %   pa = mkpp(brks,coefa);
% %   pb = mkpp(brks,coefb);
% %   pc = mkpp(brks,coefc);
% %   A0 = ppval(pa,eps2);	    
% %   B0 = ppval(pb,eps2);	   
% %   C0 = ppval(pc,eps2);	   
% % case 1,
% %  
% %   [E1, H1] = meshgrid(Tpp,Hm00);
% %   A0 = interp2(E1,H1,A00,Tp,Hm0,method);
% %   B0 = interp2(E1,H1,B00,Tp,Hm0,method);
% %   C0 = interp2(E1,H1,C00,Tp,Hm0,method);
% % end
