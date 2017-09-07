function [tc, ind, l_ind]=dat2tc(x,v,wdef)
%DAT2TC Extracts troughs and crests from data.
%
% CALL:  [TC, tc_ind, v_ind] = dat2tc(x,v,wdef);
%
%      TC     = a two-column matrix with the trough and crest turning
%                points, times in first column, values in second.
%
%      tc_ind = indices to the trough and crest turningpoints 
%               of the original sequence x.
%
%      v_ind  = indices to the level v crossings of the original 
%               sequence x. (d,u)
%
%        x    = the surface elevation.
%
%        v    = the reference level (default  v = mean of  x).
%
%	wdef   = defines the type of wave. Possible options are
%           dw', 'uw', 'tw', 'cw' or 'none'. Default is 'none'.
%           If wdef='none' all troughs and crests will be returned,
%           otherwise only the paired ones will be returned
%           according to the wavedefinition.
% 
% Example:
%   x = load('sea.dat'); x1 = x(1:200,:);
%   tc = dat2tc(x1,0,'dw');
%   plot(x1(:,1),x1(:,2),tc(:,1),tc(:,2),'ro',x1(:,1),zeros(1,200),':');
%
%   assert(length(tc), 20);
%   close all;
% 
% See also  dat2crossind,  wavedef

% Tested on: Matlab 6.0, 5.3, 5.2, 5.1
%
% History: 
% revised pab Feb2004  
% Revised by jr 02.04.2001
% - Added example, updated info. 
% Modified by Per A Brodtkorb 27.07.98
% Modified version of the old tp2tc routine which is significantly faster 
% if you use it to extract trough and crest turningpoints from the surface 
% elevation. If you use it to extract trough and crest turningpoints from a 
% sequence of turningpoints the gain in computational time is smaller. 
%
% However, this new version is more flexible. It is able to return the 
% indices to all the trough and crest turningpoints and 
% level v - crossings of the original sequence x or just those 
% turningpoints and crossings according to the wave definition. 
% (This is useful in a zero down cross analysis). 


%error(nargchk(1,3,nargin))
narginchk(1,3)
xn = x;
[n m]= size(xn);
if n<m
 b=m;m=n;n=b; 
 xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!');
end

%istime=1;

switch m
 case 1, xn=[ (1:n)' xn(:)];  %istime=0;
 case 2, % dimension OK!
 otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ');
end

if ((nargin<3) || isempty(wdef)),
  wdef='none';
end

if ((nargin<2) || isempty(v)),
  v=mean(xn(:,2));
  disp(['   The level v is set to: ', num2str(v)]);
end

%n = length(xn);

% find level v down-crossings and up-crossings
% according to wavedefinition and number of level l crossings
[l_ind Nc] = dat2crossind(xn,v,wdef);

if Nc<=2, 
  disp('Warning: There are no waves!');
  disp('         Empty matrix returned.');
  tc=[];
  return
end

% determine the number of trough2crest (or crest2trough) cycles
if mod(Nc,2),
  Ntc=(Nc-1)/2;  % Nc is odd 
else
  Ntc=(Nc-2)/2;  % Nc is even
end

% allocate variables before the loop increases the speed
tc=zeros(Nc-1,2);
ind=tc(:,1);

if (xn(l_ind(1),2)>xn(l_ind(1)+1,2)),%%%% the first is a down-crossing
  for i=1:Ntc;, 
    % trough 
    [tc(2*i-1,2) ind(2*i-1)]=min(xn((l_ind(2*i-1)+1):l_ind(2*i),2)); 
    % crest 
    [tc(2*i,2) ind(2*i)]=max(xn((l_ind(2*i)+1):l_ind(2*i+1),2)); 
  end  % for i=1:Ntc loop
  
  if (2*Ntc+1<Nc)&&(strcmp(wdef,'none')||strcmp(wdef,'tw')), 
    % trough 
    [tc(Nc-1,2) ind(Nc-1)]=min(xn((l_ind(Nc-1)+1):l_ind(Nc),2));
  end
else %%%% the first is a up-crossing
  for i=1:Ntc 
    % crest 
    [tc(2*i-1,2) ind(2*i-1)]=max(xn((l_ind(2*i-1)+1):l_ind(2*i),2)); 
    % trough 
    [tc(2*i,2) ind(2*i)]=min(xn((l_ind(2*i)+1):l_ind(2*i+1),2)); 
  end  % for i=1:Ntc loop
  
  if (2*Ntc+1<Nc)&&(strcmp(wdef,'none')||strcmp(wdef,'cw')), 
    % crest 
    [tc(Nc-1,2) ind(Nc-1)]=max(xn((l_ind(Nc-1)+1):l_ind(Nc),2));
  end
end
					

ind=l_ind(1:(Nc-1))+ind;
tc(:,1)=xn(ind,1);

%if (l_ind(Nc)+ind(Nc-1)<=n)&(wdef=='none'),
%    [tc(Nc,2) ind(Nc)]=max(xn((l_ind(Nc)+1):n,2)); % local max; probably 
%                                                    % a crest
%end






