function Snew=mkdspec(S,D,plotflag)
% MKDSPEC Make a directional spectrum
%         frequency spectrum times spreading function
%
% CALL:  Snew=mkdspec(S,D,plotflag)
%  
%       Snew = directional spectrum (spectrum struct)
%       S    = frequency spectrum (spectrum struct)
%                  (default jonswap)  
%       D    = spreading function (special struct)
%                  (default spreading([],'cos2s'))
%       plotflag = 1: plot the spectrum, else: do not plot (default 0)   
%
% Creates a directional spectrum through multiplication of a frequency
% spectrum and a spreading function: S(w,theta)=S(w)*D(w,theta)
%  
% The spreading structure must contain the following fields:
%   .S (size [np 1] or [np nf])  and  .theta (length np)  
% optional fields: .w (length nf), .note (memo) .phi (rotation-azymuth)   
%  
% NB! S.w and D.w (if any) must be identical.
%
% Example: 
%  S=jonswap
%  D=spreading(linspace(-pi,pi,51),'cos2s')
%  Snew=mkdspec(S,D,1) 
%  
% See also  spreading, rotspec, jonswap, torsethaugen  
  

% Tested on: Matlab 5.3
% History:
% revised by jr 31.03.2001 - Field added: norm. 
% by es 23.11.1999
  
if nargin<1
  Snew=demospec;
  return
end
if isempty(S)
  S=jonswap;
end
if nargin<2||isempty(D)
  D=spreading([],'cos2s');
end
if nargin<3
  plotflag=0;
end

Snew    = createspec('dir');
Snew.w  = S.w;
Snew.tr = S.tr;
Snew.h  = S.h;
Snew.phi = D.phi;
Snew.norm = S.norm;
Snew.note = S.note;
if isfield(D,'note')
  Snew.note=[Snew.note,'; ',D.note];
end
Snew.date = datestr(now);
Snew.theta = D.theta;
S.S = S.S(:);
if ~isfield(D,'w')||isempty(D.w)||length(D.w)==1
  Snew.S = D.S(:)*S.S';
else
  if length(S.w)~=length(D.w)
    error('S.w and D.w must be identical!')
  elseif any(abs(S.w-D.w)>1e-10)
    error('S.w and D.w must be identical!')
  end    
  Snew.S=D.S.*S.S(:,ones(1,length(D.theta)))';
end

if plotflag==1
  plotspec(Snew)
end
