function k=democc_markmax(proc,tp,k,diff)

%DEMOCC_MARKMAX plots load and marks a maximum
%
% CALL: k=democc_markmax(proc,tp,k,diff)
%
% Used by democc 
  
% Tested on: matlab 5.3
% History:
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version by Mats Frendahl
  
ms = 20; % markersize

n=length(tp);
if k<1
  k=1;
elseif k>n
  k=n;
end

% Check if it is a maximum, otherwise find nearest maximum.
if k~=n
  if tp(k+1,2)>tp(k,2)
    k=k+1;
  end
else
  if tp(k-1,2)>tp(k)
    k=k-1;
  end
end
  
plot(proc(:,1),proc(:,2)), hold on
plot(tp(k,1),tp(k,2),'k.','markersize',ms)
hold off

axis([min(proc(:,1)) max(proc(:,1)) 1.1*min(proc(:,2)) 1.1*max(proc(:,2))])
