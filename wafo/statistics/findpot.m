function [Ie, executionTime] = findpot(Data,Ymin,Tmin)
%FINDPOT Find indices to Peaks over threshold values
%
% CALL: [Ie] = findpot(Data,Ymin,Tmin)
%
% Ie        = indices to extreme values, i.e., all Y > Ymin which are 
%             at least Tmin distance apart.
% Te        = corresponding sampling times for Ye.
% Data      = [T, Y], a two-column matrix with sampling times and values in
%             first and second column, respectively.
% Ymin      = minimum threshold for levels in Y.
% Tmin      = minimum distance to another peak [same unit as T] (default 1)
%
% FINDPOT finds indices to peaks over threshold values, i.e., all Y > Ymin,
% which are at least Tmin distance apart.
%
% Example:
% x    = load('sea.dat'); x1 = x(1:400,:);
% tc   = dat2tc(x,0,'dw');
% ymin = 2*std(x(:,2));
% tmin = 10; % sec
% I = findpot(x,ymin,tmin);
% y = x(I,2);t = x(I,1);
% Ie = findpot(tc,ymin,tmin);
% ye = tc(Ie,2);te = tc(Ie,1);
% plot(x(:,1),x(:,2),tc(:,1),tc(:,2),'ro',x(:,1),zeros(1,length(x)),':',te, ye,'k.',t,y,'+'),shg
%
% See also fitgenpar, decluster, extremalidx


%History
% revised pab 2007
% -fixed a bug when data are too close
% by pab 2005

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


%error(nargchk(2,3,nargin))
narginchk(2,3)
if (nargin<3 || isempty(Tmin))
  Tmin = 1;
end


Ie = find(Data(:,2)>Ymin);
Ye = Data(Ie,2);
Te = Data(Ie,1);
if length(Ye)<=1
  return
end


dT = diff(Te);
notSorted = any(dT<0);
if notSorted
  [Te,I] = sort(Te);
  Ie    = Ie(I);
  Ye    = Ye(I);
  dT     = diff(Te);
end

isTooSmall = (dT<=Tmin);

if any(isTooSmall)
  tstart     = clock;

  isTooClose = [isTooSmall(1); isTooSmall(1:end-1) | isTooSmall(2:end); isTooSmall(end)];
 
  % Find opening (NO) and closing (NC) index for data beeing to close:
  iy = findextrema([0;0; isTooSmall;0]);
  
  NO = iy(1:2:end)-1;
  NC = iy(2:2:end)-1;
  
  for ix = 1:length(NO);
    iz  = NO(ix):NC(ix);
    iOK = findOKpeaks(Ye(iz),Te(iz),Tmin);
    
    if any(iOK)
%      isTooClose(iz(iOK)) = 0;
      isTooClose(NO(ix)+iOK-1) = 0;
    end
  end
  
  %% Remove data which is too close to other data.
  tooClose = find(isTooClose);
  if any(tooClose)
    Ie(tooClose) = [];
  end
  executionTime = etime(clock,tstart);
end

function iOK = findOKpeaks(Ye,Te,Tmin)
%findOKpeaks Return indices to the largest maxima that are at least Tmin
% distance apart.

Ny      = length(Ye);


[tmp,I] = sort(-Ye); % sort in descending order
%Ye1     = Ye(I);
Te1     = Te(I);
oOrder(I)  = 1:Ny; % indices to the variables original location



isTooClose = zeros(Ny,1);
  
POOL      = zeros(Ny,2);
POOL(1,:) = Te1(1) + [-Tmin,Tmin];
K         = 1;
for N = 2:Ny
  isTooClose(N) = any((POOL(1:K,1)<=Te1(N) & Te1(N)<=POOL(1:K,2)));
  if (not(isTooClose(N)))
    K         = K + 1 ;
    POOL(K,:) = Te1(N) + [-Tmin,Tmin];
  end
end
iOK = find(not(isTooClose(oOrder)));
