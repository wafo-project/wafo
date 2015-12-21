% RECLEANUP Clears global variables defined and used by the RECDEMO
%

global pwdstr RECFIGNUM

if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
if ~isempty(pwdstr)
  cd(pwdstr) % cd to the directory where recdemo was started
end  
% clear Defined Globals 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear global pwdstr recmenulabels recfilename xn map wind 
clear global inds zcrit dcrit ddcrit        % findoutliers output, input
clear global Nsim L csm1 csm2 param         % reconstruct input
clear global xr g g2 test0 tobs mu1o mu1oStd % reconstruct output                
clear global Sr m0 m2 Tm02 Hm0 Vrms Hrms    % Spectral density, moments and Sea state parameters
clear global V H rate                       % dat2steep output, input
clear global phat res noverlap              % dist2dfit output, input
clear global sphat CSMA CSMB                % dist2dsmfun2 output, input
clear global ft2 fkde                       % dist2dpdf2 and kdebin outputs
clear global kernel hs L2                   % kdebin input
RECFIGNUM=0;
%end