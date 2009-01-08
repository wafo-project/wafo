function h=trplot(g,g2,ma,sa)
%TRPLOT Plots transformation, g, eg. estimated with dat2tr.
%
% CALL:  trplot(g,g2,ma,sa)
%
%  g,g2  = smoothed and empirical estimate of the transformation g,
%  ma,sa = mean and standard deviation, respectively, of the observed
%          function. 
%          Default  ma = mean(g(:,1)) 
%                   sa = (abs(g(1,1)-ma)+abs(g(end,1)-ma))/10
%
% See also  dat2tr, lc2tr, mctrtest

% Tested on: Matlab 6.0, 5.3, 5.2, 5.1
%
% History:
% revised pab Feb2004  
% revised jr 03.04.2001
% - fixed a bug regarding nargin
% - updated information
% revised pab 04.01.2001
% - added the possibility that g is a transformation object
% revised pab 01.01.2001
% - added ih
% modified by svi 29.09.99
% g and g2 are compared with the linear transformation based on (ma,sa).
% Obs. estimates of the transformation are not normalized.
% by pab 11.11.98
%

error(nargchk(1,4,nargin))
switch class(g)
case 'double',
 if nargin<3||isempty(ma),  ma=mean(g(:,1)); end
 if nargin<4||isempty(sa),  sa=(abs(g(1,1)-ma)+abs(g(end,1)-ma))/10;end
 case 'struct' , % transformation object.
  tr = g;
  [g,ma,sa, form] = trunmak(tr); % split object
  switch form
   case 'pp', error('Not implemented for ''pp'' form yet.')
   case 'table', 
     if isempty(ma),  ma=mean(g(:,1)); end
     if isempty(sa),  sa=(abs(g(1,1)-ma)+abs(g(end,1)-ma))/10;end
  end
end

color='rgbwkymc'; 
uu=(g(:,1)-ma)/sa;
ih = ishold;
if ih,  ix=3;else  ix=1;end
  
hh = plot(g(:,1),g(:,2),color(ix),g(:,1),uu,'g--');

if nargin>1&&~isempty(g2)
  hold on
  stairs(g2(:,1),g2(:,2))  
  if ~ih, hold off, end
end

%axis([uu(1) uu(end) uu(1) uu(end)])
%axis square



title('Estimated transform')  
ylabel('g(u)')
xlabel('u')

if nargout==1,  h=hh;end

