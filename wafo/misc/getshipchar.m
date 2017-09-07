function shipchar = getshipchar(value,property)
%GETSHIPCHAR Estimate ship characteristics from value of one ship-property
%
% CALL: shipchar = getshipchar(value,property)
%
% shipchar = struct containing estimated mean values and standard-deviations
%            of ship characteristics:
%            .max_deadweight    [kkg], (weight of cargo, fuel etc.)
%            .length            [m]
%            .beam              [m]
%            .draught           [m]
%            .servicespeed      [m/s]
%            .propellerdiameter [m]
%
% value    = value (or values) to use in the estimation.
% property = string defining the ship property used in the estimation.
%           Options are:
%           'max_deadweight' (default),'length','beam','draft','servicespeed',
%           'propellerdiameter'.
%
% The deadweight is the remaining loading capacity at any time, and 
% displacement is the weight of ship and cargo at any time. Maximum
% deadweight relative to maximum displacement the is around 0.85 for
% tankers and bulk carriers, and around 0.75 for container ships. Ref: 
% ['Basic principles of ship propulsion' www.manbw.com/article_003859.html]
%
% Lightweight is the weight of an empty ship, and is typically only
% reported for military vessels. See ref. for conversion factors if
% required. (NB: NRT, GRT, NT, GT (net and gross register tons) relates to
% volume, not to tonnage. 1 GRT = 100 ft^3 = 2.83 m^3.) 
% 
% The length was found from statistics of 40 vessels of size 85 to
% 100000 tonn. An exponential curve through 0 was selected, and the
% factor and exponent that minimized the standard deviation of the relative
% error was selected. (The error returned is the same for any ship.) The
% servicespeed was found for ships above 1000 tonns only.  
%
% The propeller diameter formula is from Gray and Greeley, (1978).
%
% Example
%  ss = getshipchar(10,'servicespeed');
% 

% Other units: 1 ft = 0.3048 m and 1 knot = 0.5144 m/s


% Reference
% Gray and Greeley, (1978), 
% "Source level model for propeller blade rate radiation for the world's merchant 
% fleet", Bolt Beranek and Newman Technical Memorandum No. 458.


% by pab March 2006
% based on tonnage2geo by
% Stig Synnes, October 2005




%error(nargchk(1,2,nargin))
narginchk(1,2)
if nargin<2 || isempty(property)
  property = 'max_deadweight';
end
switch lower(property(1))
  case 'l', % length
    property = 'length';
    max_deadweight = ((value./3.45).^(2.5));
  case 'b', % beam
    property = 'beam';
    max_deadweight = ((value./1.78).^(1./0.27));
  case 'd', % draft
    property = 'draught';
    max_deadweight = ((value./0.80).^(1./0.24));
  case 's', % service speed
    property = 'servicespeed';
    max_deadweight = ((value./1.14).^(1./0.21));
  case 'p', % propeller diameter
    property = 'propellerdiameter';
    L = (value/0.12).^(4/3);
    max_deadweight = ((L./3.45).^(2.5));
  otherwise % deadweight
    property = 'max_deadweight';
    max_deadweight = value;
end
propertySTD = [property 'STD'];


L    = round(3.45.*max_deadweight.^0.40);
Lerr = L.*0.13;


W    = round(1.78.*max_deadweight.^0.27.*10)./10;
Werr = W.*0.10;

D    = round(0.80.*max_deadweight.^0.24.*10)./10;
Derr = D.*0.22;

%S    = round(2/3*(L).^0.525);
S    = round(1.14.*max_deadweight.^0.21.*10)./10;
Serr = S.*0.10;


Pdiam     = 0.12*L.^(3/4);    
Pdiam_err = 0.12*Lerr.^(3/4);

max_deadweight = round(max_deadweight);
max_deadweightSTD = 0.1*max_deadweight;

sz = size(L);

shipchar = struct('max_deadweight', max_deadweight,'max_deadweightSTD', max_deadweightSTD, ... 
  'length',L,'lengthSTD',Lerr, ...
  'beam',W,'beamSTD',Werr, ...
  'draught',D,'draughtSTD',Derr, ...
  'servicespeed',S,'servicespeedSTD',Serr, ...
  'propellerdiameter',Pdiam, 'propellerdiameterSTD',Pdiam_err );

shipchar.(propertySTD) = zeros(sz);
