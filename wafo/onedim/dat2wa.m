function [T, index, ind] = dat2wa(xn,h,pdef,wdef,rate),
%DAT2WA Extracts sequence of wavelengths from data.
%
%  CALL:  [T, index] = dat2wa(x,v/h,pdef,wdef/index,rate);
%
%      T = sequence of waveperiods (or wavelengths).	
%      x = two column data matrix with sampled times and values.
%    v/h = reference level ( default v=mean(x(:,2)) ) or
%          rainflow filtering height (default h=0)	 
%   pdef = defines the type of waveperiod (wavelength) returned:
%          Level v separated 't2c', 'c2t', 't2t' or 'c2c' -waveperiod.
%          Level v 'd2d', 'u2u', 'd2u' or 'u2d' -waveperiod.
%          Rain flow filtered (with height greater than h) 
%          'm2M', 'M2m', 'm2m' or 'M2M' -waveperiod.
%          Explanation to the abbreviations:
%          M=Max, m=min, d=down-crossing, u=up-crossing , 
%          t=trough and c=crest.
%          Thus 'd2d' means period between a down-crossing to the
%          next down-crossing and 'u2c' means period between a 
%          u-crossing to the following crest.
%          (Default is 'd2d')
%   wdef = defines the type of wave. Possible options are
%          'mw','Mw','dw', 'uw', 'tw', 'cw' or 'none'. Default is 'none'.
%          If wdef='none' all troughs and crests will be used,
%          otherwise only the troughs and crests which define a
%          wave according to the wavedefinition are used.
%
%  index = index sequence of one of the following : 
%          -level v-crossings (indices to "du" are required to 
%           calculate 'd2d', 'd2u', 'u2d' or 'u2u' waveperiods) 
%          -level v separated trough and crest turningpoints  
%           (indices to 'tc' are required to calculate 
%           't2t', 't2c', 'c2t' or 'c2c' waveperiods)
%          -level v crossings and level v separated trough and 
%           crest turningpoints (indices to "dutc" are
%           required to calculate t2u, u2c, c2d or d2t
%           waveperiods)  
%          -rainflow filtered turningpoints with minimum rfc height h 
%           (indices to "mMtc" are required to calculate
%           'm2m', 'm2M', 'M2m' or 'M2M' waveperiods)
%
%  rate = interpolation rate. If rate larger than one, then x is
%         interpolated before extrating T
%
% Example:
%  x = load('sea.dat'); x1 = x(1:400,:); 
%  [T, ind] = dat2wa(x1,0,'c2c'); % Returns crest2crest waveperiods
%  subplot(121); waveplot(x1,'-',1,1); subplot(122); histgrm(T);
%
%  assert(length(T), 22);
%  assert(ind(1:5)', [ 12, 29, 32, 40, 57]);
%  close all;
%
% See also  dat2tp, dat2tc, dat2crossind, perioddef


% This is a more flexible version than the dat2hwa or tp2wa routines. 
% There is a secret option: if pdef='all' the function returns 
% all the waveperiods 'd2t', 't2u', 'u2c' and 'c2d' in sequence.
% It is up to the user to extract the right waveperiods.
% If the first is a down-crossing then the first is a 'd2t' waveperiod.
% If the first is a up-crossing then the first is a 'u2c' waveperiod.
%
%Example:
% [T ind]=dat2wa(x,0,'all'); %returns all waveperiods
% nn = length(T);
% % want to extract all t2u waveperiods
% if x(ind(1),2)>0 % if first is down-crossing
%   Tt2u=T(2:4:nn);
% else  % first is up-crossing
%   Tt2u=T(4:4:nn);
% end
%
% assert(nn, 22);
% assert(ind(1:5), []);

% Tested on: Matlab 5.3, 5.2, 5.1
% History:
% revised pab Feb2004  
% revised pab 28.06.2001
% - added secret option ind = indices to actual points used in the
%   calculations
% - added call to ecross => improved accuracy in the level v-crossing
%   period calculations
% revised pab 03.12.1999
%  - added interpolation before extracting the parameters
% last modified by Per A. Brodtkorb 07.08.98

  error(nargchk(1,5,nargin))
x=xn;

[n m]= size(x);
if n<m
 b=m;m=n;n=b; 
 x=x';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

switch m
 case 2, % dimension OK!
 otherwise, error('Wrong dimension of input! dim must be 2xN or Nx2 ')
end

if ((nargin<4) || isempty(wdef)),
 wdef='none';
 index=[];
else
 if ischar(wdef)
   index=[];
 else
   index=wdef;
   wdef='none';
 end
end

if nargin<5||isempty(rate)||(rate<=1),  % no interpolation
  
elseif rate>1 % interpolate with spline
  xx=x;
  dT=xx(2,1)-xx(1,1);
  dT=dT/rate;
  ti=(xx(1,1):dT:xx(end,1))';
  %interpolate=1;
  x=zeros(length(ti),2);
  x(:,1)=ti;
  x(:,2)=interp1(xx(:,1),xx(:,2),ti,'*spline'); 
end

if ((nargin<3) || isempty(pdef)),
 pdef = 'd2d';
end

if ((nargin<2) || isempty(h)) ,
  if (pdef(1)=='m') || (pdef(1)=='M'),
    h=0;
    disp(['   The minimum rfc height, h,  is set to: ', num2str(h)])
  else
    h=mean(x(:,2));
    disp(['   The level l is set to: ', num2str(h)])
    % l is h in order to save number of arguments in the 
    % function call. Hopefully not too confusing!
  end
end

if ( isempty(index)),
  switch pdef(1:3) % comparing only the three first characters
    
    case {'m2m', 'm2M', 'M2m','M2M'} , % 
      %find rainflow filtered min and max turning points
      [TP index]=dat2tp(x,h,wdef);       

    case {'u2u','u2d','d2u', 'd2d'},%
      %find level v down-crossings and up-crossings
      index=dat2crossind(x,h,wdef);

    case { 't2t','t2c','c2t' 'c2c'}, %
      %find level v trough and crest turningpoints
      [TC, index]=dat2tc(x,h,wdef);
      
    case { 'd2t','t2u', 'u2c', 'c2d','all'  }, 
      %find level v down-crossings, up-crossings and trough and crest turningpoints
      [TC index l_ind]=dat2tc(x,h,wdef);
      index=(sort([index ; l_ind])); % sorting crossings and tp in sequence
    
    otherwise, error('Unknown option!');
  end 
end % if nargin

if (x(index(1),2)>x(index(2),2)), % if first is down-crossing or max 
  switch pdef
    case 'all', start=1;	% secret option!
    case {'d2t','M2m','c2t', 'd2u' , 'M2M','c2c','d2d'}, start=1;
    case {'t2u','m2M', 't2c', 'u2d' ,'m2m','t2t','u2u'}, start=2;
    case 'u2c', start=3;
    case 'c2d', start=4;
    otherwise, error('Unknown option!');
  end %switch pdef
else 			% first is up-crossing or min 
  switch pdef 
    case 'all', start=1;	% secret option!
    case {'u2c','m2M', 't2c', 'u2d','m2m','t2t','u2u'} , start=1;
    case {'c2d','M2m','c2t', 'd2u' , 'M2M','c2c','d2d'}, start=2;
    case 'd2t', start=3;
    case 't2u', start=4;
    otherwise, error('Unknown option!');
  end % switch pdef
end %if 

% determine the steps between wanted periods
switch pdef
  case {'d2t','t2u','u2c', 'c2d' }, step=4; 
  case 'all', step=1;% secret option!  
  otherwise, step=2; 
end % switch pdef

% determine the distance between min2min, t2t etc..
switch pdef
  case {'m2m','t2t','u2u','M2M','c2c','d2d'}, dist=2;
  case 'all', dist=1;% secret option!
  otherwise,  dist=1;
end % switch pdef

nn = length(index);
% New call: (pab 28.06.2001)
if strmatch(pdef(1),{'u','d'}),
  t0 = ecross(x(:,1),x(:,2),index(start:step:(nn-dist)),h);
else % min, Max, trough, crest or all crossings wanted
  t0 = x(index(start:step:(nn-dist)),1); 
end
if strmatch(pdef(3),{'u','d'}),
  t1 = ecross(x(:,1),x(:,2),index((start+dist):step:nn),h);
else % min, Max, trough, crest or all crossings wanted
  t1 = x(index((start+dist):step:nn),1); 
end
T = t1-t0;
if nargout>2, % Secret option: indices to the actual crossings used.
  index=index(:);
  ind = [index(start:step:(nn-dist)) index((start+dist):step:nn)].';
  ind = ind(:);
end

return

% Old call: kept just in case
%T  = x(index((start+dist):step:nn),1)-x(index(start:step:(nn-dist)),1);












