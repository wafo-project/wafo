function  [RFC,RFC1,res,def] = tp2rfc(x,def,RFC0,res0)
%TP2RFC Finds the rainflow cycles from the sequence of turning points.
%
% CALL:  [RFC,RFC1,res] = tp2rfc(tp,def,RFC0,res0);
%                   RFC = tp2rfc(tp);
%
% Output:
%   RFC     = Rainflow cycles (residual included).      [N,2]/[N,4]
%   RFC1    = Rainflow cycles (without resudual).       [N1,2]/[N1,4]
%   res     = Residual.                               [nres,1]/[nres,2]
%
% Input:
%   tp      = Turning points.                            [T,1]/[T,2]
%   def     = Choice of definition of rainflow cycles   [struct array]
%   def.res = Treatment of residual.
%             'up':   Count min-to-Max cycles,    (default)
%                     gives correct number of upcrossings.
%             'down': Count Max-to-min cycles, 
%                     gives correct number of downcrossings.
%             'CS':   Cloormann/Seeger method, 
%                     gives all closed hysteresis loops.
%                     This method is identical to the French AFNOR recommendation, 
%                     and the ASTM standard (variant called simplified version).
%   def.time = 0: Don't store time of max and min. (default)
%              1: Store the time when the maxima and minima occured in columns 3-4. 
%   def.asymmetric = 0: gives the symmetric RFC (default),
%                    1: gives the asymmetric RFC (or From-To RFC), time order between 
%                       maximum and rainflow minimum is preserved.
%   RFC0    = Rainflow cycles (without resudual).       [N0,2]/[N0,4]
%   res0    = Residual.                               [nres0,1]/[nres0,2]
%
% Calculates the rainflow cycles (RFC) for the sequence of turning points, 
% by using the so-called 4-point algorithm.
%
% It is possible to split the signal into smaller parts, and calculate 
% RFC part by part. It can be especially useful for long signals.
% We count the first part and for the second part we continue counting 
% from previously counted 'RFC0' with residual 'res0':
%   [RFC1,RFC0,res0] = tp2rfc(tp(1:1000,:));      % Firts 1000 points
%   [RFC2] = tp2rfc(tp(1001:end,:),[],RFC0,res0); % Point 1001 to end
% This shall give the same result as (i.e. ARFC=ARFC2)
%   [RFC] = tp2rfc(tp);                           % Calculate all at once
%   sum(RFC~=RFC2)                                % Shall return  [0 0]
%
% This routine doesn't use MEX, Fortran or C codes, only matlab code.
%
% Example:
%   x = load('sea.dat'); tp=dat2tp(x);
%   RFC1 = tp2rfc(tp);      % Default (min-to-Max cycles in residual)
%   ccplot(RFC1); 
%   RFC2 = tp2rfc(tp,'CS'); % Compare with AFNOR/ASTM standard
%   [I,J] = find(RFC1(:,1)~=RFC2(:,1) | RFC1(:,2)~=RFC2(:,2));
%   hold on; plot(RFC1(I,1),RFC1(I,2),'b+',RFC2(I,1),RFC2(I,2),'rx'); hold off;
%
%   close all;
%
% See also  findrfc, dat2tp, rfcfilter, tp2arfc
  
% Further examples:
%   % Rainflow cycles with time
%   def.res='up'; def.time=1; % Store times 
%   RFC=tp2rfc(tp,def); RFC(1:10,:), ccplot(RFC)
%
%   % For long signals it is possible to split the input in smaller parts
%   [dummy,RFC0,res0] = tp2rfc(dat2tp(x(1:5000,:)));     % First part
%   [RFC3] = tp2rfc(dat2tp(x(5001:end,:)),[],RFC0,res0); % Second part
%   % RFC3 shall be the same as RFC1. Check this!
%   ccplot(RFC1), hold on,plot(RFC3(:,1),RFC3(:,2),'r.'), hold off 

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (P�r Johannesson) 2000-01-04
%   Uses same syntax as 'tp2rfc' in WAT
% Revised by PJ 26-Jul-2000
%   New format of def.
%   Added input 'RFC0' and 'res0'. New output 'RFC1' and 'res'
%   Now supports AFNOR and ASTM standards for rainflow counting.
% Revised by PJ 06-Jul-2005
%   Fixed error with def & mod to avoid warning i R14SP2.

% Check input arguments
ni = nargin;
no = nargout;
%error(nargchk(1,4,ni));
narginchk(1,4)
if ni < 2, def=[]; end
if ni < 3, RFC0=[]; end
if ni < 4, res0=[]; end

% Check input def
def0=def;
if ~isempty(def)
  if isnumeric(def)
    def=[]; def.time = def0;
  elseif ischar(def)
    def=[]; def.res = def0;
  elseif ~isstruct(def)
    def=[];
  end
end

% Set default values
if ~isfield(def,'res')
  def.res = 'up';
end
if ~isfield(def,'time')
  def.time = 0;
end
if ~isfield(def,'asymmetric')
  def.asymmetric = 0;
end

% Count rainflow cycles
if no<2
  ARFC = tp2arfc(x,def,[],res0);
else
  [ARFC,ARFC1,res] = tp2arfc(x,def,[],res0);
end

% Convert to symmetric RFC ?
if def.asymmetric == 0 % Symmetric rainflow cycles
  RFC = make_symmetric(ARFC);
else
  RFC = ARFC;
end

% Add previously counted cycles (if any)
if ~isempty(RFC0)
  RFC = [RFC0; RFC];
end

% Rainflow cycles without residual
if no>2, 
  if def.asymmetric == 0 % Symmetric rainflow cycles
    RFC1 = make_symmetric(ARFC1);
  else
    RFC1 = ARFC1;
  end
  % Add previously counted cycles (if any)
  if ~isempty(RFC0)
    RFC1 = [RFC0; RFC1];
  end
end


function RFC = make_symmetric(ARFC)

  I = ARFC(:,1)>ARFC(:,2);
  [N,M]=size(ARFC);
  if M == 2 % No time
    J=1;
  else      % Time of occurances is stored in column 3:4
    J=[1 3];
  end
  RFC = ARFC;
  RFC(I,J) = ARFC(I,J+1);
  RFC(I,J+1) = ARFC(I,J);
