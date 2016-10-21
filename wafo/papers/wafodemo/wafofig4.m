function wafofig4
% WAFOFIG4  Directional spectra using Thorsethaugen and cos2s spreading 
%           function with s=15 (solid) and s frequency dependent (dash)
%

  
global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end
global St Tp ma mb sp
D=spreading(51,'cos2s',0,[sp sp 2*pi/Tp 0 0],St.w,0);
Std=mkdspec(St,D); 
Dw=spreading(linspace(-pi,pi,51),'cos2s',0,[sp sp 2*pi/Tp ma mb],St.w);
Stdw=mkdspec(St,Dw); 

plotspec(Std,1,'k-'), hold on
plotspec(Stdw,1,'k--'),hold off
wafostamp('Figure 4','(ER)')    
