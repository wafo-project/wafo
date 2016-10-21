function recfig3
% RECFIG3 Example of the reconstructed (solid) and original data set (pluses)
%         Use the mouse to zoom in on a specific region
%

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global xn xr m0 inds
%choose a random place to see how it looks
Ns=ceil(rand*length(xr));
%if (rand<27/39),
% Ns=ceil(rand*26000);
%else
% Ns=ceil(30000+rand*8700);
%end

waveplot(xr,'k-',xn(inds,:), 1,1,sqrt(m0),4)
zoom xon
%axis([xr(Ns,1)/60 xr(Ns+100,1)/60 -7 7])
%waveplot(xr(20201:20300,:),xn(inds,:), 1,1,sqrt(m0),4)
wafostamp('Figure 3','(NR)')
