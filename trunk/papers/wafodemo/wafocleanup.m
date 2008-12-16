% WAFOCLEANUP Clears global variables defined and used by the WAFODEMO
%

global pwdstr WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  return
end
if ~isempty(pwdstr)
  cd(pwdstr) % cd to the directory where wafodemo was started
end  
% clear Defined Globals 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear global pwdstr wafomenulabels xn Jxn Nxr map wind 
clear global Fs                          % sampling frequency for torsethaugen
clear global St Hm0 Tp                   % torsethaugen output, input Fig 2
clear global xt Nsim Iseed               % spec2sdat output , input    Fig1
clear global Ste L                       % dat2spec output, input      Fig2
clear global fTt fTte u Np nit speed     % spec2thpdf output, input    Fig3  
clear global Tt rate                     % dat2wa output, input        Fig3
clear global ma mb sp                    % spreading input             Fig4

clear global fTcfAc NNp Nh Nnit Nspeed   % spec2thpdf output, input    Fig5  
clear global NVcf NHd Nrate              % dat2steep output, input     Fig5

clear global fTcfAcTc JNp Jh Jnit Jspeed % spec2thpdf output, input    Fig6  
clear global JTcf JAc Jind Jrate         % dat2steep output, input     Fig6

clear global kdeTt kernel hs L2          % kdebin input: kernel, smoothing  (Fig3)
clear global kdeVcfHd Nkernel Nhs NL2    % and transformation  parameter,   (Fig5)
clear global kdeTcfAcTc Jkernel Jhs JL2  % respectively                     (Fig6)

WAFOFIGNUM=0;
%end