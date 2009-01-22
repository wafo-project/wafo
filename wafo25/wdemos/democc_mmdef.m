function ccmM=democc_mmdef(proc,tp,point,ccpt)

%DEMOCC_MMDEF illustrates the definition of min-max cycles.
%
% CALL: ccmM=democc_mmdef(proc,tp,point,ccpt)
%
% Used by democc 

% Tested on: matlab 5.3
% History:
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version by Mats Frendahl

ms = 20; % markersize

democc_markmax(proc,tp,point,0);
hold on, title('Definition of peak-trough count')
plot(tp(point-1,1),tp(point-1,2),'k.','markersize',ms)
plot(tp(point,1),tp(point,2),'k.','markersize',ms)
hold off
if isempty(ccpt)
   ccmM=[tp(point-1,2) tp(point,2)];
else
   ccmM=[ccpt; tp(point-1,2) tp(point,2)];
end


