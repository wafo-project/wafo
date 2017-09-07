function [indx, indy] = wfindpeaks(S,Np,min_h,min_p)
%FINDPEAKS find peaks of vector or matrix possibly rainflow filtered
%  
% CALL: [ix iy] = wfindpeaks(S,Np,min_h,min_p)   
%
%    ix     = row subscripts to the peaks of S if iy present
%             otherwise linear index to peaks of S
%    iy     = column subscripts to the peaks of S
%    S      = matrix or vector
%    Np     = The Np highest peaks are found (if exist). (default 2)
%    min_h  = The threshold in the rainflowfilter (default 0.05*range(S(:))).
%             A zero value will return all the peaks of S.
%    min_p  = 0..1, Only the peaks that are higher than 
%                   min_p*max(max(S))  min_p*(the largest peak in S)
%                   are returned (default  0).
% Example: 
%    x = (0:0.01:10); 
%    S = x.^2+10*sin(3*x)+0.5*sin(50*x); 
%    clf; plot(x,S);
% % Find highest 8 peaks that are not less that 0.3*"global max" and have 
% % rainflow amplitude larger than 5.
%    assert(wfindpeaks(S',8,5,0.3), [909, 695, 482]');
%
%   close all;
% 
% See also  dat2tp

% GL, Feb 2011, Changed name to  wfindpeaks  to avoid name collision 
 % Tested on: matlab 5.3, 8, 9
% History:
% revised pab Feb2004
% -changed default value for min_h
% -added nargchk  
% revised pab 12.10.1999
%   fixed a bug  
% by Per A. Brodtkorb  25.09.1999
  
%error(nargchk(1,4,nargin));
narginchk(1,4)
if (nargin <4) || isempty(min_p)
  min_p=0;
end
if (nargin <3) || isempty(min_h)
  min_h = 0.05*range(S(:));
end

if (nargin <2) || isempty(Np)
 Np = 2;
end
Ssiz=size(S);
if min(Ssiz)<2
  ndim=1;
  S=S(:)';
  m=max(Ssiz);
  n=1;
else 
  ndim=2;
  n=Ssiz(1);
  m=Ssiz(2);
end

x=(1:m).';

% Finding turningpoints of the spectrum
% Returning only those with rainflowcycle heights greater than h_min 
indP=[]; % indices to peaks
ind=[];
for iy=1:n % find all peaks
  TuP =dat2tp([x S(iy,:)'],min_h); % find turning points
  
  if ~isempty(TuP)
    ind=TuP(2:2:end,1); % extract indices to maxima only
  else % did not find any , try maximum
    [y ind]=max(S(iy,:));
  end
  if ndim>1,
    if iy==1, ind2 = find((S(iy,ind)>S(iy+1,ind)) );
    elseif iy==n, ind2 = find((S(iy,ind)>S(iy-1,ind)));
    else,
      ind2 = find((S(iy,ind)>S(iy-1,ind)) & (S(iy,ind)>S(iy+1,ind)));
    end
    if ~isempty(ind2)
       indP = [indP;ind(ind2) iy*ones(length(ind2),1)];
    end      
  end
end

if isempty(ind)||(isempty(indP)&&(ndim>1))
  indx=[];
  indy=[];
  return
end
if ndim>1
  % indP(1) = col subscript, indP(2) = row subscript
  ind=sub2ind([n m],indP(:,2),indP(:,1));
end

[Y ind2] = sort(S(ind));
ind2 = flipud(ind2(:)) ; 

% keeping only the Np most significant peak frequencies.
ind=ind(ind2(1:min(Np,length(ind))));
if (min_p >0 ), % Keeping only peaks larger than min_p percent relative
                % to the maximum peak 
  [Y ind2]=find(S(ind) >min_p*max(S(ind)));
  ind=ind(ind2);
end

if nargout>1
  [indx indy]=ind2sub(Ssiz,ind);
else
  indx=ind;
end
