function cocc(param,cc,matrix,clevels,fs)
%COCC Plots cycles as points together with isolines of a cycle matrix.
%
% CALL: cocc(param,cc,matrix)
%       cocc(param,cc,matrix,clevels,psize)
%
%        param   = the parameter matrix.
%        cc      = a two column matrix with cycles.
%        matrix  = an  nxn  matrix.
%        clevels = a vector with levels of isolines.
%        psize   = point size, (defult value is 12).
%
% Plots cycles of min- and max-values as a point process in the plane
%  together with isolines of an input matrix.
%
% Example:
%   x=load('sea.dat');
%   tp = dat2tp(x);
%   rfc = tp2rfc(tp,'CS');           % Rainflow cycles
%   RFM = dat2rfm(x,0,[-3 3 100]);  % Rainflow matrix
%   ERFM = rfmextrapolate(RFM);     % Extrapoalted RFM
%   cocc([-3 3 100],rfc,RFM);
%   cocc([-3 3 100],rfc,ERFM);
%
% See also  cmatplot, plotcc


% Copyright (c) 2003 by Pär Johannesson

% Tested  on Matlab  6.5
%
% History:
% Created from old version of WAFO (Jan-2001)
% Updated by PJ 03-Jun-2003

% Check input arguments

if nargin<5
  fs=[];
end

if nargin<4  
  clevels=[];
end

if isempty(fs)
  fs=12;
end

if isempty(clevels)
  fmax=max(max(matrix));
  clevels=fmax*[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8];
end

F = matrix';
u=levels(param);
contour(u,u,F,clevels,'r')
axis('square')  

hold on
if ~isempty(cc)
  plot(cc(:,1),cc(:,2),'.','markersize',fs)  % PJ 970415
end
plot(param(1:2),param(1:2),'k--')
hold off
xlabel('min')
ylabel('Max')

 clevels=sort(clevels);
 n_clevels=length(clevels);
 if n_clevels>12
   disp('   Only the first 12 levels will be listed in table.')
   n_clevels=12;
 end

textstart_x=0.65;
textstart_y=0.45;
delta_y=1/33;
h=figtext(textstart_x,textstart_y,'Level curves at:','normalized');
set(h,'FontWeight','Bold')

textstart_y=textstart_y-delta_y;

for i=1:n_clevels
  textstart_y=textstart_y-delta_y;
  figtext(textstart_x,textstart_y,num2str(clevels(i)),'normalized')
end

