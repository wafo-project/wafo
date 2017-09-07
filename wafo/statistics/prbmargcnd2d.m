function  [y ,eps1] = prbmargcnd2d(phat,x1lo,x1up,x2lo,x2up)
%PRBMARGCND2D returns the probability for rectangular regions.
%
% CALL: [P tol] = prbmargcnd2d(phat,x1lo,x1up,x2lo,x2up);
%
%   P    = probability
%   tol  = absolute tolerance, i.e., abs(int-intold)
%   phat = parameter structure (see fitmargcnd2d)
%   xilo = lower integration limits
%   xiup = upper integration limits
% 
%  The size of P is the common size of XILO and XIUP.  
% 
% Example
%  x1=linspace(0,10)';
%  phat.x = {[x1,exp(-0.1*x1)] 2 };
%  phat.dist = {'rayl','rayl'};
%  prb = prbmargcnd2d(phat,1,2,1,2);
%  f = pdfmargcnd2d2(x1,x1,phat);
%  pdfplot(f); hold on;
%  plot([ 1 1 2 2 1],[1 2 2 1 1]); hold off;
%
%  See also  fitmargcnd2d rndmargcnd2d pdfmargcnd2d cdfmargcnd2d

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



% tested on: matlab 5.2
% history:
% revised pab 27.10.2000
%  - added example text
%  Per A. Brodtkorb 28.10.98

%error(nargchk(5,5,nargin))
narginchk(5,5)
%defining global variables
global PHAT CONDON
condon=CONDON; % save old value
CONDON=1;
if (nargin < 5), 
  error('Requires 5 input arguments.'); 
end
eps2=1e-5;%relative tolerance
% nit toolbox function
[y eps1] = gaussq2d('cdfmargcnd2dfun',x1lo,x1up,x2lo,x2up,eps2);
CONDON=condon; %restore the default value


