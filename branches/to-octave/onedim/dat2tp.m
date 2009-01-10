function [tp, ind] = dat2tp(x,h,wdef)
%DAT2TP Extracts turning points from data,
%        optionally rainflowfiltered. 
%
% CALL:  [tp ind] = dat2tp(x,h,wdef);
%
%    x  = two column data matrix with sampled times and values.
%
%    tp = a two column matrix with times and turning points.
%
%   ind = indices to the turning points in the original sequence.
%
%    h  = a threshold; 
%         if  h<0, then  tp=x; 
%         if  h=0, then  tp  is a sequence of turning points (default); 
%         if  h>0, then all rainflow cycles with height smaller than
%                  h  are removed.
%
%  wdef = defines the type of wave. Possible options are
%	  'mw' 'Mw' or 'none'. (Default 'none').
%         If wdef='none' all rainflow filtered min and max 
%         will be returned, otherwise only the rainflow filtered 
%         min and max  which define a wave according to the 
%         wave definition will be returned.
%
% Example:
%   x  = load('sea.dat'); x1 = x(1:200,:);
%   tp = dat2tp(x1,0,'Mw'); tph = dat2tp(x1,0.3,'Mw');
%   plot(x1(:,1),x1(:,2),tp(:,1),tp(:,2),'ro',tph(:,1),tph(:,2),'k*')
% 
%  See also  findcross, findrfc, tp2rfc, dat2tc

% Tested on: matlab 6.0, 5.3, 5.2, 5.1

% History:
% revised pab 2008
% -wdef option now works correctly
% revised pab oct 2005
% -replaced some code with a call to findextrema
% revised pab Feb2004  
% Revised by jr 03.04.2001
% - added example, updated info 
% Modified by Per A. Brodtkorb 07.08.98
% This is a modified version which is about 20 to 30 times faster than  
% the version of dat2tp in WAT (performance on a pentiumII 233 MHz  
% with 32 MB ram and Matlab 5.0 under Linux). The reason is
% that this version does not save x to disk. Instead it passes 
% the arguments directly to the executeable file. 
% This new version is also more flexible. It is able to return the 
% indices to the turningpoints 
% (This is useful when determining the wave steepness etc...). 

error(nargchk(1,3,nargin))

xn = x;
[n m]= size(xn);
if n<m
  b=m;m=n;n=b; 
  xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

%istime=1;

switch m
  case 1, x2=xn; %istime=0;
  case 2, x2=xn(:,2);% dimension OK!
  otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ')          
end

if ((nargin<3) || isempty(wdef)),
  wdef='none';
end

if (nargin<2) || isempty(h),
   h=0;
end

if h<0 
  tp=xn; 
  ind=(1:n)';
  disp('Warning: h<0, the original data is returned')
  return   
end			

ind = findextrema(x2);

if length(ind)<2, 
  tp=[];
  return;
end

% In order to get the exact up-crossing intensity from rfc by 
% mm2lc(tp2mm(rfc))  we have to add the indices
% to the last value (and also the first if the 
% sequence of turning points does not start with a minimum).

if  x2(ind(1))>x2(ind(2)),		
 % adds indices to  first and last value
 ind=[1; ind ;n]; 
else %adds index to the last value
 ind=[ind; n];
end

if h>0 
   ind1 = findrfc(x2(ind),h);
   ind  = ind(ind1);
end 


Nm=length(ind); % number of min and Max


switch wdef  % switch wdef
 case {'mw','Mw'},
   % make sure that the first is a Max if wdef == 'Mw'
   % or make sure that the first is a min if wdef == 'mw'
   if xor((x2(ind(1))>(x2(ind(2)))),strcmp(wdef,'Mw')),
     ind=ind(2:Nm);
     Nm=Nm-1;
   end

   % make sure the number of minima and Maxima are according to the wavedef.
   % i.e., make sure Nm=length(ind) is odd
   if ~(mod(Nm,2)), % if Nm is even do
     ind(Nm)=[];
     Nm=Nm-1;
   end	
   
 case {'none'}% do nothing
 otherwise, error('Unknown wave definition') 
end

tp=xn(ind,:);

