function democc_plotmat(demow2,ccrfc,ccmM)
%DEMOCC_PLOTMAT plots RFC and min-max counts.
%
% CALL: democc_plotmat(demow2,ccrfc,ccmM)
%
% Used by democc  
  
% Tested on: matlab 5.3
% History:
% Created by by PJ 10-Jan-2000
  
figure(demow2)

subplot(1,2,1)
if ~isempty(ccrfc)
  plotcc(ccrfc),hold on
  v=axis; plot(v(1:2),v(3:4),'--'),hold off
end
title('Rainflow cycles')

subplot(1,2,2)
if ~isempty(ccmM)
  plotcc(ccmM),hold on
  v=axis; plot(v(1:2),v(3:4),'--'),hold off
end
title('min-max cycles')
