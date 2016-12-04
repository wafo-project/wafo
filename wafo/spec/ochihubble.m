function [S, S1] = ochihubble(w,sdata,plotflag)
%OCHIHUBBLE Calculates (and plots) a bimodal Ochi-Hubble spectral density
% 
%  CALL:  [S S1] = ochihubble(w,data,plotflag); 
%
%        S    = a struct containing the total spectral density
%        S1   = a 2D struct array containing the swell and wind generated 
%               part of the spectral density
%        w    = angular frequency (default linspace(0,3,257))
%        data = [Hm0 def]
%               Hm0 = significant wave height (default 7 (m))
%               def = defines the parametrization (default 1)
%                     1 : The most probable spectrum
%                     2,3,...11 : gives 95% Confidence spectra
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the  spectrum.
%
%  The OH spectrum is a six parameter spectrum, all functions of Hm0.
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
% Example
% S = ochihubble;
% plotspec(S)
%
% See also  mkochihubble, jonswap, torsethaugen

% References:
% Ochi, M.K. and Hubble, E.N. (1976)
% 'On six-parameter wave spectra.'
% In proc. 15th Conf. Coastal Engng.,1,pp301-328

% Tested on: matlab 5.2
% History:
% revised pab aug 2007
% - changed name from ohspec2 to ochihubble
% - replaced call to ohspec with a call  to mkochihubble.
% by pab 16.02.2000

monitor=0;

if nargin<3||isempty(plotflag)
  plotflag=0;
end
if nargin<1||isempty(w)
  w=linspace(0,3,257).';
end
sdata2=[7 1]; % default values
if nargin<2||isempty(sdata)
else
 ns1=length(sdata);
 k=find(~isnan(sdata(1:min(ns1,2))));
 if any(k)
   sdata2(k)=sdata(k);
 end
end %

Hm0=sdata2(1);
def=sdata2(2);
 

[H,Hs,Hw] = mkochihubble('Hm0',Hm0,'def',def);
S      = createspec;
S.S = H(w);
S.w = w(:);
S.norm = 0;
S.note = sprintf('Ochi-Hubble, Hm0 = %5.2g, def = %d',Hm0,def);
if nargout>1 || plotflag>0
  S1      = createspec;
  S1(1).S = Hs(w);
  S1(2).S = Hw(w);
  [S1.w]  = deal(w(:));
  [S1.norm] = deal(0); % The spectrum is not normalized
  opt1 = Hs('options');
  opt2 = Hw('options');
  S1(1).note = sprintf('Ochi-Hubble, Hm0 = %5.2g, Tp = %5.2g',opt1.Hm0,opt1.Tp);
  S1(2).note = sprintf('Ochi-Hubble, Hm0 = %5.2g, Tp = %5.2g',opt2.Hm0,opt2.Tp);
end



if monitor
  disp(['Hm0, def      = ' num2str([Hm0 def])]);
end


if plotflag
  ih=ishold;
  plotspec(S,plotflag);
  hold on
  plotspec(S1,plotflag,'k--');
  if ~ih,hold off,end
end
