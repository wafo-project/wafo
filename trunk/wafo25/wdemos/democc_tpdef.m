function democc_tpdef(L)

%DEMOCC_TPDEF illustrates the definition of turning points.
%
% CALL: democc_tpdef(L)
%
% Used by democc 

% Tested on: matlab 5.3
% History:
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version by Mats Frendahl

k=1; hold on
plt = plot(0,min(L(:,1)),'.','erasemode','xor','markersize',16);
pltmax = plot(L(k,1),L(k,2),'c.','erase','none','markersize',10);
pltmin = plot(L(k,1),L(k,2),'r.','erase','none','markersize',10);
for k = 2:length(L)-1
   set(plt,'xdata',L(k,1),'ydata',L(k,2))
   if ( (L(k-1,2)<=L(k,2)) && (L(k,2)>L(k+1,2)) )
      set(pltmax,'xdata',L(k,1),'ydata',L(k,2))
   elseif ( (L(k-1,2)>=L(k,2)) & (L(k,2)<L(k+1,2)) )
      set(pltmin,'xdata',L(k,1),'ydata',L(k,2))
   end
   drawnow
end;
hold off
