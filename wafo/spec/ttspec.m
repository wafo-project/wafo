function Snew = ttspec(S,varargin)
%TTSPEC Toggle Transform between angular frequency and frequency spectrum
%
% CALL:  Snew = ttspec(S,ftype,thtype)
%
%     Snew = spectrum with new frequency type
%        S = spectrum 
%    ftype = 'f' if Snew should be given with frequency in Hz
%            'w' if Snew should be given with angular frequency in rad/s
%            (default is the opposite of what is given in S)
%   thtype = 'radians' if Snew.theta should be given in radians (default)
%            'degrees' if Snew.theta should be given in degrees 
%            
%
% If ftype is the same frequency type as for S then Snew = S
% If S is a wavenumber spectrum then Snew = S.
% 
% Note: The order of ftype and thtype is arbitrary and only the
%       letters 'f','w','r' or 'd' is needed for unique identification
%
% Examples: % Change from angular frequency to frequency in Hz and from
%           % angle in radians to angle in degrees.  
%    S   = demospec('dir');
%    Sf  = ttspec(S);
%    Sf1 = ttspec(Sf,'f','d'); % = ttspec(Sf,'f','degrees'); 
%    plotspec(S,3), figure(2)
%    plotspec(Sf,3), figure(3)
%    plotspec(Sf1,3)
%
%    close('all')
%
% See also  datastructures



%Tested on: Matlab 5.3, 5.2
% History: 
% revised pab 22.06.2001
% - fixed a bug: S.phi field was not checked when 'd', or 'r' option was used.
% revised pab 13.06.2000
%  - added more checks on input
%  - added thtype
%  - removed recursive call to itself => more efficient code.
% revised pab 08.02.2000
%  - S can now be an array of structs
% revised pab 24.01.2000
%  -added ftype
% by pab 13.11.99



% Error checking
%~~~~~~~~~~~~~~~
if ~isstruct(S)
  error('Input must be a spectral density struct)')
end
if ~(strcmpi(S(1).type(end-2:end),'req') || strcmpi(S(1).type(end-2:end),'dir'))
  disp('This is not a frequency spectrum. Nothing is changed.')
  Snew = S; % return old if not a freq
  return
end

ind = [isfield(S,'f'), isfield(S,'w')];
if all(ind==0),
  error('This is not a correct spectral density struct: w and f field does not exist')
elseif all(ind==1),
  error('This is not a correct spectral density struct: w and f field can not both exist')
end

% Setting default values
%~~~~~~~~~~~~~~~~~~~~~~~~~
ftypeold = 'wf';
ftype    = ftypeold(ind); 
ftypeold = ftypeold(~ind);
thtype   = 'radians'; 

P = varargin; Np = length(P);
for ix=1:Np,
  switch lower(P{ix}(1)),
    case {'f','w'}, ftype = lower(P{ix});
    case {'d','r'}, thtype = lower(P{ix});
  end
end


% Main computations
%~~~~~~~~~~~~~~~~~~~
isphi = isfield(S , 'theta');
Ns       = length(S(:));
if strcmp(ftypeold,ftype)
  Snew = S; % nothing is changed
  
  if isfield(S , 'theta')
    for iy=1:Ns,
      if strcmpi(S(iy).type(end-2:end),'dir'),  % Directional spectrum 
	if (abs(max(S(iy).theta)-min(S(iy).theta))>3*pi)  % theta given in degrees
	  if strcmpi(thtype(1),'r'),  % radians wanted
	    Snew(iy).theta = S(iy).theta*pi/180;
	    Snew(iy).S     = S(iy).S*180/pi;
	    if isphi,  Snew(iy).phi     = S(iy).phi*pi/180;end
	  end 
	else                             % theta given in radians
	  if strcmpi(thtype(1),'d'),  % degrees wanted 
	    Snew(iy).theta = S(iy).theta*180/pi; 
	    Snew(iy).S     = S(iy).S*pi/180;
	    if isphi,  Snew(iy).phi     = S(iy).phi*180/pi;end
	  end 
	end
      end
    end
  end
else % change to new type
  fnames   = fieldnames(S(1));
  ind      = strmatch(ftypeold,fnames);
  fnames   = {fnames{[1:ind-1, ind+1:end] }};
  
  Snew  = createspec(S(1).type,ftype);
  for ix = 1:length(fnames),             % make sure all non-standard elements of S
                                         % are also transfeered to Snew
    Snew.(fnames{ix}) = [];
  end
  Snew(Ns) = Snew(1); % make sure it is an array of structs
  
  for iy=1:Ns,
    for ix = 1:length(fnames),
      Snew(iy).(fnames{ix}) = S(iy).(fnames{ix});
    end
    if strcmp(ftypeold,'f')    
      Snew(iy).w = S(iy).f*2*pi;
      Snew(iy).S = S(iy).S/(2*pi);
    else
      Snew(iy).f = S(iy).w/(2*pi);
      Snew(iy).S = S(iy).S*(2*pi);
    end
    
    if strcmpi(S(iy).type(end-2:end),'dir'),  % Directional spectrum 
      if (abs(max(S(iy).theta)-min(S(iy).theta))>3*pi)  % theta given in degrees
	if strcmpi(thtype(1),'r'), % radians wanted
	  Snew(iy).theta = S(iy).theta*pi/180; 
	  Snew(iy).S     = Snew(iy).S*180/pi;
	   if isphi,  Snew(iy).phi     = S(iy).phi*pi/180;end
	end 
      else                            % theta given in radians
	if strcmpi(thtype(1),'d'), % degrees wanted 
	  Snew(iy).theta = S(iy).theta*180/pi;
	  Snew(iy).S     = Snew(iy).S*pi/180;
	   if isphi,  Snew(iy).phi     = S(iy).phi*180/pi;end
	end 
      end
    end
  end % iy
end

