function [S, S1] = ohspec3(w,def,plotflag)
% OHSPEC3 Calculates Bimodal Ochi-Hubble spectral densities
%         from mixed sea states
%
%  CALL:  [S S1] = ohspec3(w,def,plotflag); 
%
%        S    = a struct containing the total spectral density
%        S1   = a 2D struct array containing the swell and wind generated 
%               part of the spectral density
%        w    = angular frequency (default linspace(0,3,257))
%        def  = vector or a scalar containing numbers between 1 and 9. 
%               1,2,3 : Swell dominated Sea state
%               4,5,6 : Wind sea dominated  Sea state
%               7,8,9 : Mixed wind-sea and swell systems with comparable energy.
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the  spectrum.
%
%  The OH spectrum is a six parameter spectrum.
%  Here 9 different parameterizations are defined representing
%  3 different type of sea states. Each of the of the 3 sea states are 
%  parameterized with 3 different inter-modal distances. The exact parameters
%  are given below where the subscripts w = wind sea,  s = swell and
%  Hm0 = significant wave height, Tp = peak period,  L = spectral shape parameter
%
%             Target spectra parameters:
%    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   def     Hm0s    Hm0w     Tps     Tpw     Ls     Lw
%  -----------------------------------------------------
%    1      5.5     3.5     14.3     9.1     3.0    6.5
%    2      6.5     2.0     14.3     6.7     3.5    4.0
%    3      5.5     3.5     22.2     6.5     3.0    6.0
%    4      2.0     6.5     14.3     9.1     3.0    6.5
%    5      2.0     6.5     14.3     6.7     4.0    3.5
%    6      2.0     6.5     22.2     6.5     2.0    7.0
%    7      4.1     5.0     14.3     9.1     2.1    2.5
%    8      4.1     5.0     14.3     6.7     2.1    2.5
%    9      4.1     5.0     22.2     6.5     2.1    2.5 
%
% Example: % Create an array struct with all the parameterizations
%    
%     S = ohspec3([],1:9);
%  
% See also  bretschneider, jonswap, torsethaugen

% References:
% Ochi, M.K. and Hubble, E.N. (1976)
% 'On six-parameter wave spectra.'
% In proc. 15th Conf. Coastal Engng.,1,pp301-328
%
% Rodriguez, G.R. and Soares, C.G. (2000)
% "Wave period distribution in mixed sea states"
% In proc. ETCE/OMAE2000 Joint Conference Energy for the new millenium.
%
% Rodriguez, G.R., Soares, C.G. and Pacheco, M. (2004)
% "Wave period distribution in mixed sea states"
% J. Offshore mech. Arct. Eng. Vol 126, pp 105-112

% Tested on: matlab 5.2
% History:
% revised pab Apr2005
% fixed a bug: ohspec([],1:9,1) works correctly
% by pab 16.02.2000

monitor=0;

if nargin<3||isempty(plotflag)
  plotflag=0;
end
if nargin<1||isempty(w)
  w=linspace(0,3,257).';
end
if nargin<2||isempty(def)
 def =1;
elseif length(def)>1
  for ix=1:length(def),
    [S(ix) S1(ix,1:2)]=ohspec3(w,def(ix));
  end %
  if plotflag
     plotspec(S,plotflag);
  end
  return
end

switch def
  % Swell dominated sea states
  case 1, Hm0i=[5.5 3.5]; Tpi=1./[0.07 0.11];Li=[3   6.5];
  case 2, Hm0i=[6.5 2.0]; Tpi=1./[0.07 0.15];Li=[3.5 4];
  case 3, Hm0i=[5.5 3.5]; Tpi=1./[0.045 0.155];Li=[3 6];
    % Wind dominated sea states
  case 4, Hm0i=[2 6.5]; Tpi=1./[0.07 0.11];Li=[3 6.5];
  case 5, Hm0i=[2 6.5]; Tpi=1./[0.07 0.15];Li=[4 3.5];
  case 6, Hm0i=[2 6.5]; Tpi=1./[0.045 0.155];Li=[2 7];
    % Mixed Wind sea and swell systems with comparable energy
  case 7, Hm0i=[4.1 5]; Tpi=1./[0.07 0.11];Li=[2.1 2.5];
  case 8, Hm0i=[4.1 5]; Tpi=1./[0.07 0.15];Li=[2.1 2.5];
  case 9, Hm0i=[4.1 5]; Tpi=1./[0.045 0.155];Li=[2.1 2.5]; 
  otherwise
    Hm0i=[5.5 3.5]; Tpi=1./[0.07 0.11];Li=[3   6.5];
    def=1;
end

if 0, % the following make Hm0 the same for all the cases
  
  switch def
    % Swell dominated sea states
  case 1, Hm0i=[5.5 3.5]*0.99184053150473; Tpi=1./[0.07 0.11];Li=[3   6.5];
  case 2, Hm0i=[6.5 2.0]*0.95078090199753; Tpi=1./[0.07 0.15];Li=[3.5 4];
  case 3, Hm0i=[5.5 3.5]*0.99184046072521; Tpi=1./[0.045 0.155];Li=[3 6];
    % Wind dominated sea states
  case 4, Hm0i=[2 6.5]*0.95078087556934; Tpi=1./[0.07 0.11];Li=[3 6.5];
  case 5, Hm0i=[2 6.5]*0.95078087556934; Tpi=1./[0.07 0.15];Li=[4 3.5];
  case 6, Hm0i=[2 6.5]*0.95078087556934; Tpi=1./[0.045 0.155];Li=[2 7];
  end
end

Ni   = 4*Li+1;
Mi   = [4 4];

Hs = mkbretschneider('Hm0',Hm0i(1),'Tp',Tpi(1),'N',Ni(1),'M', Mi(1));
Hw = mkbretschneider('Hm0',Hm0i(2),'Tp',Tpi(2),'N',Ni(2),'M', Mi(2));

S      = createspec;
Hs1 = Hs(w);
Hw1 = Hw(w);
S.S = Hs1+Hw1;
S.w = w(:);
S.norm = 0;
Hm0 = 4*sqrt(trapz(w,S.S));
S.note = sprintf('Ochi-Hubble, Hm0 = %5.2g, def = %d',Hm0,def);
if nargout>1 || plotflag>0
  S1      = createspec;
  S1(1).S = Hs1;
  S1(2).S = Hw1;
  [S1.w]  = deal(w(:));
  [S1.norm] = deal(0); % The spectrum is not normalized
  opt1 = Hs('options');
  opt2 = Hw('options');
  S1(1).note = sprintf('Ochi-Hubble, Hm0 = %5.2g, Tp = %5.2g',opt1.Hm0,opt1.Tp);
  S1(2).note = sprintf('Ochi-Hubble, Hm0 = %5.2g, Tp = %5.2g',opt2.Hm0,opt2.Tp);
end


% for ix=1:2,
%   S1(ix)=ohspec(w,[Hm0i(ix), Tpi(ix),Li(ix)]);
% end
% S=S1(1);
% S.S=S1(1).S+S1(2).S;
% 
% %S.note=['Ochi-Hubble3, Hm0 = ' num2str(Hm0i)  ', Tp = ' num2str(Tpi) , ' L = ' num2str(Li)];
% S.note=['Ochi-Hubble3, def = ' num2str(def) ];
% if monitor
%   disp(['Hm0, Tp      = ' num2str([Hm0i Tpi])])
% end

if plotflag
  ih=ishold;
  plotspec(S,plotflag);
  hold on
  plotspec(S1,plotflag,'k--');
  if ~ih,hold off,end
end
