%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% testhspdfs.m --- 
%% Author          : Per Andreas Brodtkorb
%% Created On      : Tue Feb 17 19:17:02 2004
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Tue Feb 17 19:26:03 2004
%% Update Count    : 7
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sc = 0.25;
hc = 3;

% Conditional probability of steep and high waves given seastates
% i.e., Prob(Hd>hc,Scf>sc|Hs,Tz)  
upperTail = 1;
Hs = linspace(2.5,11.5,30);
Tp = linspace(4.5,19.5,48);
[T,H] = meshgrid(Tp,Hs); 

pnl = thsnlcdf(hc,sc,H,T,upperTail);
pl  = thscdf(hc,sc,H,T,upperTail);
p87 = mk87cdf(hc,sc,H,T,upperTail);
 
global THSNLPAR
Tp0  = THSNLPAR.Tp;
Hm00 = THSNLPAR.Hm0;
Tm020 = THSNLPAR.Tm02;
[Tp1,Hs1] = meshgrid(Tp0,Hm00);
Tm02Nl = interp2(Tp1,Hs1,Tm020,T,H,'cubic');


v = 10.^(-6:-1);  
contourf(Tp,Hs,log10(pnl),log10(v))
xlabel('Tp')
ylabel('Hs')  
fcolorbar(log10(v))  