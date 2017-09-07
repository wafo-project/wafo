function [I,R] = sudg(n,p,varargin)
%SUDG Some Useful Design Generators
%
% CALL:  [I, R] = sudg(n,p);
%
% I = matrix of generating relations. size p X q
% R = Resolusion
% n = number of variables.
% p = number of generators.
%
% SUDG may be used in conjunction with FFD to construct two-level
% fractional factorial designs with the highest possible resolution.
% In general, a 2^(n-p) fractional design is produced by P generators and
% has a defining relation containing 2^P words. The resolution R of a
% fractional design is the length of the shortest word in the defining
% relation. A resolution R design has a complete factorial (possibly
% replicated) in every subset of R-1 variables. In general, a design of
% resolution R is one in which no k-factor effect is confounded with any
% other effect containing less then R-k factors.
%
% See also ffd, cl2cnr, cnr2cl, yates, fitmodel, getmodel, alias, cdr, plotresponse, nplot

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




% Reference 
% Box, G.E.P, Hunter, W.G. and Hunter, J.S. (1978)
% Statistics for experimenters, John Wiley & Sons, chapter 12,pp 410


% Tested on: Matlab 5.3
% History:
% By Per A. Brodtkorb 16.03.2001

%error(nargchk(2,3,nargin))
narginchk(2,3)
nmp = n-p;
if isempty(p)||p<=0, I = zeros(0,1); R = nmp+1; return;end


if p ==1,
  I = 1:n;
  R = nmp+1;
elseif n+1==2^nmp, % Saturated Resolution III design
  I = zeros(n,nmp);
  iz = 0;
  for ix = 1:nmp,
    iz      = iz+1;
    I(iz,1) = ix;
    iz0     = iz;
    for iy = 1:iz0-1,
      iz = iz+1;
      I(iz,:)   = I(iy,:);
      ind        = find(I(iy,:)==0, 1 );
      I(iz,ind) = ix;
    end
  end
  
  k = find(I(:,2)~=0); % Find only interactions
  I = [I(k,:) (nmp+1:n).'];
  R = 3;
else
  % The following is taken from Box et-al(1978) pp 410
 txt =  'Requested design is not currently available or impossible';
  switch nmp,
    case 3, %8,
      if p>4, error(txt), end      
      I0 =  [1:2 0; 1 3 0; 2:3 0; 1:3];
      R =3;
    case 4, %16,
      if p>11, error(txt), end
      I0 = [1:3 0; 2:4 0; 1 3:4 0; 1 2 4 0; 1:4; 1:2 0 0; ....
	    1 3 0 0;1 4 0 0; 2 3 0 0; 2 4 0 0; 3 4 0 0 ];
      if p>=5, R=3;else R=4;end
    case 5, %32,
      if p>6, error(txt), end
      if p== 6;
	I0 = [1:3;2:4;3:5;1 3:4; 1 4:5; 2 4:5];
      elseif p==3
	I0 =[1:3 0; 1:2 4 0; 2 3:5];
      else
	I0 = [ 1 2 3 5 ; 1 2 4 5; 1 3 4 5; 2 3  4 5;  1 2 3 4];
      end
      R = 4;
    case 6,% 64,
      if p>5, error(txt), end
      if p==5,
	I0 = [1:4; 3:5 0; 2 4 5:6; 1 4:6; 1:2 6 0];
      elseif p==4
	I0 = [2:4 6; 1 3:4 6; 1 2 4:5; 1 2 3 5 ];
      else
	I0 = [1:4; 1:2 5:6; 2 4:6  ];
      end
      if p>=3, R=4;else R=5;end
    case 7, %128,
      if p>4, error(txt), end
      if p==2,
	I0 = [1 3:4 6:7; 2:3 5:7 ];
	R = 6;
      else
	I0 = [1:3 7 0 0 0; 2:5 0 0 0; 1 3:4 6 0 0 0; 1:7 ];
	R = 5;
      end
      
    otherwise,
      error(txt)
  end
  I = [I0(1:p,:) (nmp+1:n).'];
end
I = sort(I,2);
if nargin<3 && n<=50; % secret option
  I = cnr2cl(I); % Convert column number to column label
end




