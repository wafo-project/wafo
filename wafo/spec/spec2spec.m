function [Snew]=spec2spec(S,newtype,phi,v)
%SPEC2SPEC Transforms between different types of spectra
%
% CALL:  Snew = spec2spec(S,newtype,phi,v)
%
%     Snew = spectrum with new type
%        S = spectrum 
%  newtype = new type given as a string, the allowed types: 
%            'dir','freq','k2d','k1d','encdir','enc',
%            (default S.type, i.e. no change)
%          
%  In special cases, specify also:
%      phi = rotation angle of the cordinate system (x,y) counter-clockwise 
%            (rotation of spectrum clockwise), e.g. direction of ship
%             when encountered spectrum is to be computed; the wavenumber
%             spectrum along a line going through origin  and with
%             angle phi (default 0).
%        v = speed of ship, when transforming to an encountered spectrum
%            (default 0)
%  
% SPEC2SPEC performs a change of variables from the type given in input
% spectrum S, to the given 'newtype', when possible.
% Remark1: the field .note is not changed, but may need an update by hand.
% Remark2: the encountered spectrum is derived by means of numerical
% variable change, hence dense spectrum is needed to get good accuracy of
% the transformation. For negative velocities the enc. spectrum has a
% singularity hence spec2mom may give uncorrect values for the moments.
%  
% Example: 
%    S    = demospec('dir');
%    Snew = spec2spec(S,'enc',pi/6,10);
%
% See also  datastructures, rotspec

% Tested on: Matlab 5.3  
% History: 
% revised by pab 21.09.2004:
%  -Moved code to rotspec + small cosmetics fixes   
% revised by IR 27.06.2004: implemented freq -> enc.
% revised by IR 03.04.2001: Major changes.
% revised by jr 31.03.2001
% - case 'freq' -> 'enc': Changed text in error message
% revised by es 05.06.00: encdir -> freq results in enc
% revised by es 25.05.00: return if nargin<2 (default for newtype)  
% revised by es 24.05.00: removed Message: New type identical to old type
% revised by es 09.02.2000: more corrections of change to 'k1d'
% revised by es 28.01.2000: correction and improvement of change to 'k1d'  
% by es 13.08.99

%error(nargchk(1,4,nargin));
narginchk(1,4)
if nargin<0||isempty(S)
  error('Needs an input spectrum');
end

Snew = S; 
if (nargin<2||isempty(newtype))
  return
end
if (nargin<3 || isempty(phi))
  phi=0.;
end
if (nargin<4||isempty(v))
  v=0;
end

Snew = rotspec(Snew,phi);

newtype=lower(newtype);
if strcmpi(S.type,newtype)
  %   disp(' Message: New type identical to old type')
  return
end

if ~any(strcmpi(newtype,{'dir','freq','k2d','k1d','encdir','enc'}))    
  % Check if newtype is a proper type
  error('Not known new type, check spelling');
end

indim=sum(size(S.S)>1); %dimension of old spectrum
if indim==1 %then new spectrum can not have dimension larger
  if strcmp(newtype(end-2:end),'dir')||strcmp(newtype(end-2:end),'k2d')
    error('Impossible to transform from one dimension to two');
  end
end



switch lower(S.type)
  case 'dir'
   switch newtype
    case 'freq' % from 'dir'
     Snew.S = simpson(Snew.theta(:),Snew.S); 
     %integrate out angle but if the spectrum is
     %true directional spectrum then
     %one should just take SS.S(1,:)
        Snew=rmfield(Snew,'theta');
        Snew.type='freq';
      case 'k2d' % from 'dir'
        Snew=time2spa(Snew);
      case 'k1d' % from 'dir'
        Snew=spec2spec(Snew,'k2d');
        Snew=spec2spec(Snew,'k1d');
      case 'encdir' % from 'dir'
        Snew=dir2enc(Snew,2,v);
      case 'enc' % from 'dir'
        Snew=dir2enc(Snew,1,v);
    end   
   case 'freq'
    switch newtype
      case 'k1d' % from 'freq'
        Snew=time2spa(Snew);
      case 'enc' % from freq
          Snew=dir2enc(Snew,1,v);
      otherwise
        error('Specified transformation not possible');
    end
  case 'k2d'
   switch newtype
    case 'dir' % from 'k2d'
     Snew=spa2time(Snew);
    case 'freq' % from 'k2d'
     Snew=spec2spec(Snew,'dir'); %transform via dir (k1d better???)
     Snew=spec2spec(Snew,'freq');
     
    case 'k1d' 
     % from 'k2d'
     % For spectrum in Cartesian representation:
     %    The grid is rotated, the size of it is preserved
     %    (maybe it must be increased such that no nonzero points are
     %    affected, but this is not inplemented yet: i.e. corners are cut off)
     % The spectrum is assumed to be zero outside original grid.

     % PAB 2004: old calls kept just in case:
     %[k,k2] = meshgrid(S.k,S.k2);
     %[th,r] = cart2pol(k,k2);
     %phi=Snew.phi;
     %[k,k2] = pol2cart(th+phi,r);
     %Sn     = interp2(S.k,S.k2,S.S,k,k2);
     %Sn(isnan(Sn))=0.;
     %Snew.S = Sn;
     
     physicallyRotateGrid = 1;
     Snew = rotspec(Snew,0,physicallyRotateGrid);
     
     Snew.S = simpson(Snew.k2(:),Snew.S).'; %integrate out second variable
     Snew   = rmfield(Snew,'k2');
     Snew.type = 'k1d';
     if (Snew.k(1)<0)
       Snew.S = Snew.S(Snew.k>=0)+flipud(Snew.S(Snew.k<=0));
       Snew.S = Snew.S';
       Snew.k = Snew.k(Snew.k>=0);
     end
     Snew.type='k1d';
    case 'encdir' % from 'k2d'
     Snew=spec2spec(Snew,'dir');
     Snew=spec2spec(Snew,'encdir',v);
    case 'enc' % from 'k2d'
     Snew=spec2spec(Snew,'dir');
     Snew=spec2spec(Snew,'enc',v);
    otherwise 
     error('Specified transformation not possible');
   end
 case 'k1d'
  switch newtype
   case 'freq' % from 'k1d'
    Snew = spa2time(Snew);
   case 'enc' % from 'k1d'
    Snew = spa2time(Snew);
    Snew = spec2spec(Snew,'enc',v);
   otherwise      
    error('Specified transformation not possible');
  end

 case 'encdir'
  switch newtype
   case 'enc' % from 'encdir'
    Snew.S = simpson(Snew.theta(:),Snew.S); % integrate out direction
    Snew.type = 'enc';
    Snew      = rmfield(Snew,'theta');
   otherwise
    error('Specified transformation not possible or not yet available');
  end
 case 'enc'
  switch newtype 
   case 'freq' % from 'enc'  
    Snew.type='freq';
   otherwise
    error('Specified transformation not possible or not yet available');
  end
end
Snew.date=datestr(now);
