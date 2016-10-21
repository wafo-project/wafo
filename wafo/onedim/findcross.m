function ind = findcross(x,v)
%FINDCROSS Finds indices to level v up and downcrossings of a vector
%
%  CALL: ind = findcross(x,v);
%
%        x  = vector with sampled values.
%
%        v  = level v. (Default 0). 
%
%	ind = indices to the crossings in the original sequence x.
%
% Example
%  v = 0.75
%  t = linspace(0,7*pi,250); x = sin(t);
%  ind = findcross(x,v)
%  plot(t,x,'.',t(ind),x(ind),'r.', t, ones(size(t))*v)
%
% See also dat2crossind, crossdef

% there is also a mex version of this which is much faster,
% which is run instead if compiled and located before this file 
% in the MATLAB search path.

% Tested on: Matlab 5.3, 5.2 5.1

% History:
% revised pab April2004  
% fixed a bug: Now gives correct answer when turning points on crossing level
% revised pab Feb2004
% revised pab 13.06.2001
% -Added example
% - fixed a bug: this .m file previosly only returned zero crossings.
% by pab 17.07.1997


if nargin<2||isempty(v),
  v = 0;
end

xn  = int8(sign(x(:)-v));
ind = [];
n  = length(xn);
if n>1
  
  iz = find(xn==0);
  if any(iz)
    % Trick to avoid turning points on the crossinglevel.
    if iz(1)==1
      if length(iz) ==n
        warning('WAFO:FINDCROSS','All values are equal to crossing level!')
        return
      end
      
      diz = diff(iz);
      ix  = iz(find(diz>1,1,'first'));
      if ~any(ix)
        ix = iz(end);
      end
      %x(ix) is a up crossing if  x(1:ix) = v and x(ix+1) > v. 
      %x(ix) is a downcrossing if x(1:ix) = v and x(ix+1) < v. 
      xn(1:ix) = -xn(ix+1); % 
      iz(1:ix) = [];
    end


    for ix=iz.',
      xn(ix) = xn(ix-1);
    end
  end

  % indices to local level crossings ( without turningpoints)
  ind = find(xn(1:(n-1)).*xn(2:n) < 0);
 
end
end


