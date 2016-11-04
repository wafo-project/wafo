function [H,Hw,Hs] = mkochihubble(varargin)
%MKOCHIHUBBLE Function handle to bimodal spectral density.
% 
% CALL: [H,Hw,Hs] = mkochihubble(options,par1,val1,par2,val2,...)
%
%  H,Hw,Hs = function handles to to spectrum for total, wind 
%            and swell, respectively. 
%  options = structure with the fields:
%   .Hm0   = significant wave height (default 7 (m))
%   .def   = integer defining the parametrization (default 1)
%            1 : The most probable spectrum
%            2,3,...11 : gives 95% Confidence spectra
%
% MKOCHIHUBBLE Return function handle to a bimodal spectrum modelled as  
% S(w) = Ss(w) + Sw(w) where Ss and Sw are modified Bretschneider 
% spectra for swell and wind peak, respectively.
%
% The OH spectrum is a six parameter spectrum, all functions of Hm0.
%  The values of these parameters are determined from a analysis of data
%  obtained in the North Atlantic. The source of the data is the same as
%  that for the development of the Pierson-Moskowitz spectrum, but
%  analysis is carried out on over 800 spectra including those in
%  partially developed seas and those having a bimodal shape. From a
%  statistical analysis of the data, a family of wave spectra consisting
%  of 11 members is generated for a desired sea severity (Hm0) with the
%  coefficient of 0.95. 
%  A significant advantage of using a family of spectra for design  of
%  marine systems is that one of the family members yields the largest
%  response such as motions or wave induced forces for a specified sea
%  severity, while another yields the smallest response with confidence
%  coefficient of 0.95. 
%
% Example: 
%  S = mkochihubble('def',2); 
%  fplot(S,[0, 5])
%
%  close()
%
% See also  mkbretschneider, mkjonswap, mktorsethaugen
  
% References:
% Ochi, M.K. and Hubble, E.N. (1976)
% 'On six-parameter wave spectra.'
% In Proc. 15th Conf. Coastal Engng., Vol.1, pp301-328

% Tested on: matlab 6.0, 5.3
% History:
% By pab jan 2007
% based on old ohspec2 function


error(nargchk(0,inf,nargin))
options = struct('Hm0',7,'def',1,'chkseastate','on');
if (nargin==1) && ischar(varargin{1}) && strcmpi(varargin{1},'defaults')
  H = options;
  return
end
options = parseoptions(options,varargin{:});
if 1<options.def && options.def < 12
else
  options.def = 1;
end

if strcmpi(options.chkseastate,'on')
  chkseastate(options)
end

[Hw,Hs] = mkohspec(options);

H = @(w)ochihubble(options, Hw, Hs, w);
return
end %mkochihubble

%% Nested function
function S = ochihubble(options, Hw, Hs, w)
  %OCHIHUBBLE spectral density
  
  if strncmpi(w,'options',1)
    S = options;
    S.Hwoptions = Hw(w);
    S.Hsoptions = Hs(w);
    
  else 
    if options.Hm0>0
      S = Hw(w)+ Hs(w);
    else
      S = zeros(size(w));
    end
  end
end %ochihubble


%% Subfunctions

function [Hwind,Hswell] = mkohspec(options)
Hm0 = options.Hm0;
switch options.def, % 95% Confidence spectra
  case 2, hp=[0.84 0.54];wa=[.93 1.5]; wb=[0.056 0.046];Li=[3.00 2.77*exp(-0.112*Hm0)];
  case 3, hp=[0.84 0.54];wa=[.41 .88]; wb=[0.016 0.026];Li=[2.55 1.82*exp(-0.089*Hm0)];
  case 4, hp=[0.84 0.54];wa=[.74 1.3]; wb=[0.052 0.039];Li=[2.65 3.90*exp(-0.085*Hm0)];
  case 5, hp=[0.84 0.54];wa=[.62 1.03];wb=[0.039 0.030];Li=[2.60 0.53*exp(-0.069*Hm0)];
  case 6, hp=[0.95 0.31];wa=[.70 1.50];wb=[0.046 0.046];Li=[1.35 2.48*exp(-0.102*Hm0)];
  case 7, hp=[0.65 0.76];wa=[.61 0.94];wb=[0.039 0.036];Li=[4.95 2.48*exp(-0.102*Hm0)];  
  case 8, hp=[0.90 0.44];wa=[.81 1.60];wb=[0.052 0.033];Li=[1.80 2.95*exp(-0.105*Hm0)];  
  case 9, hp=[0.77 0.64];wa=[0.54 .61];wb=[0.039 0.000];Li=[4.50 1.95*exp(-0.082*Hm0)];  
  case 10,hp=[0.73 0.68];wa=[.70 0.99];wb=[0.046 0.039];Li=[6.40 1.78*exp(-0.069*Hm0)];    
  case 11,hp=[0.92 0.39];wa=[.70 1.37];wb=[0.046 0.039];Li=[0.70 1.78*exp(-0.069*Hm0)];     
  otherwise % The most probable spectrum
    %def=1;
    hp=[0.84 0.54]; wa= [.7 1.15];wb=[0.046 0.039];  Li=[3 1.54*exp(-0.062*Hm0)];
end
Hm0i = hp*Hm0;
Tpi  = 2*pi.*exp(wb*Hm0)./wa;
Ni   = 4*Li+1;
Mi   = [4 4];

Hswell = mkbretschneider('Hm0',Hm0i(1),'Tp',Tpi(1),'N',Ni(1),'M', Mi(1));
Hwind  = mkbretschneider('Hm0',Hm0i(2),'Tp',Tpi(2),'N',Ni(2),'M', Mi(2));

end % mkohspec

function chkseastate(options)

Hm0 = options.Hm0;

if Hm0<0
  error('WAFO:MKOCHIHUBBLE','Hm0 can not be negative!')
end

if options.def<1 || 11<options.def
  error('WAFO:MKOCHIHUBBLE','Def must be an integer freom 1 to 11!')
end


if Hm0==0
  warning('WAFO:MKOCHIHUBBLE','Hm0 is zero!')
end

end % chkseastate


