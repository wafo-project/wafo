function H=mkbretschneider(varargin)
%MKBRETSCHNEIDER Function handle to BRETSCHNEIDER spectral density
%
% CALL:  H = mkbretschneider(options,par1,val1,par2,val2,...)
%
%   H      = function handle to BRETSCHNEIDER spectral density. 
%  options = structure with the fields:
%   .Hm0   = significant wave height (default 7 (m))
%   .Tp    = peak period             (default 11 (sec))
%   .N      = scalar defining decay of high frequency part.   (default 5)
%   .M      = scalar defining spectral width around the peak. (default 4)
%
% MKBRETSCHNEIDER Return function handle to a BRETSCHNEIDER spectrum defined as
%
%         S(w) = A * G0 * wn^(-N)*exp(-N/(M*wn^M))
%    where 
%         G0  = Normalizing factor related to Bretschneider form
%         A   = (Hm0/4)^2 / wp     (Normalization factor)
%         wn  = w/wp       
%         wp  = 2*pi/Tp,  angular peak frequency
%
%  The BRETSCHNEIDER spectrum is a suitable model for fully developed sea, 
%  i.e. a sea state where the wind has been blowing long enough over a 
%  sufficiently open stretch of water, so that the high-frequency waves have
%  reached an equilibrium. In the part of the spectrum where the frequency is
%  greater than the peak frequency (w>wp), the energy distribution is
%  proportional to w^-5.
%  The spectrum is identical with ITTC (International Towing Tank 
%  Conference), ISSC (International Ship and Offshore Structures Congress) 
%  and Pierson-Moskowitz, wave spectrum given Hm0 and Tm01. It is also identical
%  with JONSWAP when the peakedness factor, gamma, is one.
%  For this spectrum, the following relations exist between the mean
%  period Tm01 = 2*pi*m0/m1, the peak period Tp and the mean
%  zero-upcrossing period Tz:
%
%           Tm01 = 1.086*Tz, Tp = 1.408*Tz and Tp=1.2965*Tm01
%
% Example:
%  S = mkbretschneider('Hm0',6.5,'Tp' ,10); 
%  fplot(S,[0,4])
%  options = S('options'); % get options used
%  assert(fieldnames(options), {'Hm0', 'Tp', 'N', 'M', 'chkseastate'}')
%  assert(struct2cell(options), {6.5,10,5,4,'on'}')
%
%  close()
%
% See also  mkjonswap, mktorsethaugen

% References:
%

% Tested on: matlab 6.0, 5.3
% History:
% By pab jan 2007
% New enhanced implementation based on old bretchneider function.
% renamed from bretschneider to mkbretschneider since it now returns a
% function handle instead of the old spectral struct.
% -replaced code with call to ggamspec
% -change name from pmspec to bretschneider because Bretschneider presented
% his spectrum in 1959, ITTC, ISSC and Pierson-Moskovitz spectra came first
% in 1964.

%error(nargchk(0,inf,nargin));
narginchk(0,inf)
options = struct('Hm0',7,'Tp',11,'N',5,'M',4,'chkseastate','on');
if (nargin==1) && ischar(varargin{1}) && strcmpi(varargin{1},'defaults')
  H = options;
  return
end
options = parseoptions(options,varargin{:});

if strcmpi(options.chkseastate,'on')
  chkseastate(options);
end


H = @(w)bretschneider(options, w);
return
end %mkbretschneider
%% Subfunctions

function S = bretschneider(options, w)
%BRETSCHNEIDER spectral density
%
% CALL        S = bretschneider(w)
%       options = bretschneider('options')
%
  if strncmpi(w,'options',1)
    S = options;
  else
    if options.Hm0>0
      wp = 2*pi/options.Tp;
      wn = w./wp;
      S = (options.Hm0/4)^2/wp * ggamspec(wn,options.N,options.M);
    else
      S = zeros(size(w));
    end
  end
end %bretschneider


function chkseastate(options)
  Tp = options.Tp;
  Hm0 = options.Hm0;

  if Hm0<0
    error('WAFO:MKBRETSCHNEIDER','Hm0 can not be negative!');
  end

  if Tp<=0
    error('WAFO:MKBRETSCHNEIDER','Tp must be positve!');
  end


  if Hm0==0
    warning('WAFO:MKBRETSCHNEIDER','Hm0 is zero!');
  end

end % chkseastate
