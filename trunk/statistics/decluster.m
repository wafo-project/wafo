function [Ye,Te] = decluster(Data,Ymin,Tmin)
%DECLUSTER Decluster peaks over threshold values
%
% CALL: [Ye,Te] = decluster(Data,Ymin,Tmin)
%
% Ye        = extreme values, i.e., all Y > Ymin which are 
%             at least Tmin distance apart.
% Te        = corresponding sampling times for Ye.
% Data      = [T, Y], a two-column matrix with sampling times and values in
%             first and second column, respectively.
% Ymin      = minimum threshold for levels in Y.
% Tmin      = minimum distance to another peak [same unit as T] (default 1)
%
% DECLUSTER extracts peaks over threshold values, i.e., all Y > Ymin,
% which are at least Tmin distance apart.
%
% Example:
% x    = load('sea.dat'); x1 = x(1:400,:);
% tc   = dat2tc(x1,0,'dw');
% ymin = 2*std(x(:,2));
% tmin = 10; % sec
% [ye, te] = decluster(tc,ymin,tmin);
% plot(x1(:,1),x1(:,2),tc(:,1),tc(:,2),'ro',x1(:,1),zeros(1,length(x1)),':',te,ye,'k.')
%
% See also fitgenpar, findpot, extremalidx


%History
% revised pab 2007
% -fixed a bug when data are too close
% by pab 2005
error(nargchk(2,3,nargin))
if (nargin<3 || isempty(Tmin))
  Tmin = 1;
end
Ie = findpot(Data,Ymin,Tmin);

Ye = Data(Ie,2);
Te = Data(Ie,1);
