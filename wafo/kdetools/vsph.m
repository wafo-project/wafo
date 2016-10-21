function c=vsph(d,r0)
% VSPH calculates volume of  d-dimensional sphere with radius r0
%
% CALL: vd = vsph(d,r0)
%
%  d  = dimension of sphere
%  r0 = radius of sphere (default 1)
% 

% Reference
%  Wand,M.P. and Jones, M.C. (1995)
% 'Kernel smoothing'
%  Chapman and Hall, pp 105

% History   
% revised pab dec2003
% -removed some code  
% By pab 199?  
  
  
if nargin<2 || isempty(r0)
 r0=1; % default unit sphere
end

c=(r0^d)* 2*pi^(d/2)/(d*gamma(d/2)); % Wand and Jones pp 105
return
