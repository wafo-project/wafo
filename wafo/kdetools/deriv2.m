function k=deriv2(x,y,d)
%DERIV2  High order partial derivatives of the Gaussian kernel.
%
% CALL:  k = deriv2(x,y,dstr)
%
%  k   =  partial derivatives of the 2D Gaussian kernel  
%         at the point (X,Y).
% x,y  = evaluation points
% dstr = string defining the type of partial derivative
%
% Example: 4'th p. derivative wrt. x and 2'nd p. derivative wrt. y at (0,0)
%
%    k42 = deriv2(0,0,'42');
%    assert(k42, -0.477464829275686)
%
% See also  deriv

%tested on: matlab 5.3
%revised pab 16.10.1999
%  updated to matlab 5.x + documentation
% from kdetools by   Christian C. Beardah 1995



xd=str2double(d(1));
yd=str2double(d(2));

switch xd,
  case 0,
    xterm=1;
  case 2,
    xterm=x.^2-1;
  case 4,
    xterm=x.^4-6*x.^2+3;
  case 6,
    xterm=x.^6-15*x.^4+45*x.^2-15;
  case 8,
    xterm=x.^8-28*x.^6+210*x.^4-420*x.^2+105;
  otherwise
   error('Unsupported order of derivative')     
end;
switch yd,
  case 0,
    yterm=1;
  case 2,
    yterm=y.^2-1;
  case 4,
    yterm=y.^4-6*y.^2+3;
  case 6,
    yterm=y.^6-15*y.^4+45*y.^2-15;
  case 8,
    yterm=y.^8-28*y.^6+210*y.^4-420*y.^2+105;
  otherwise
    error('Unsupported order of derivative')
end;

k=xterm.*yterm.*mkernel(x,y,'gauss');