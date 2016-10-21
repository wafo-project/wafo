%DEMOCC A program for visualization of cycle counts in random
%   loads. 
%
%   In Demonstration Window 1, the realisation is shown.
%   It is possible to mark the turning points (TP) and, for local maxima
%   chosen by the user, find rainflow cycles and min-max
%   cycles. In Demonstration Window 2, the cycle counts are
%   illustrated; peak-trough cycles to the left, rainflow cycles to
%   the right.
%  
%   NB! A realisation of a random process called 'proc' must exist in
%   workspace. 
%
% Example:
%   x = load('sea.dat');
%   proc=x(1:200,:);
%   democc

% Tested on: matlab 5.3
% History:
% Original version by Mats Frendahl
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Revised by JR 10-July-2000
%   line 24-27. Matrix dimensions
% Revised by PJ 13-Jun-2003
%   Change figure positions
% Updated by PJ 07-Jul-2005
%   Added example

democc_demow2=figure('Name','The rainflow & peak-trough cycle count','NumberTitle','off','Position',[0 300 500 350]);
democc_demow1=figure('Name','Demonstration Window 1','NumberTitle','off','Position',[0 0 1000 300]);
clf; 
democc_time=1:length(proc); 
if min(size(proc))==2, % Size of simulated proc may be 
  proc=proc(:,2);      % nx1, 1xn or nx2. 
end 
democc_y=[democc_time(:) proc(:)];         
democc_L=democc_y(:,2); democc_n=length(democc_L); clc, clf, subplot(1,1,1)

democc_ccrfc=[]; democc_ccmM=[]; democc_tp=dat2tp(democc_y); 

set(gca,'box','on'), xlabel('time'), ylabel('load')
democc_k=democc_markmax(democc_y,democc_tp,1,0);
democc_F = uicontrol('style','push','units','normal','pos',[.92 .93 .03 .06], ...
        'string','+1','call','democc_k=democc_k+2; democc_k=democc_markmax(democc_y,democc_tp,democc_k,-2);');
democc_FF= uicontrol('style','push','units','normal','pos',[.96 .93 .03 .06], ...
        'string','+5','call','democc_k=democc_k+10; democc_k=democc_markmax(democc_y,democc_tp,democc_k,-10);');
democc_REW = uicontrol('style','push','units','normal','pos',[.92 .86 .03 .06], ...
        'string','-1','call','democc_k=democc_k-2; democc_k=democc_markmax(democc_y,democc_tp,democc_k,2);');
democc_REW = uicontrol('style','push','units','normal','pos',[.96 .86 .03 .06], ...
        'string','-5','call','democc_k=democc_k-10; democc_k=democc_markmax(democc_y,democc_tp,democc_k,10);');
democc_RFC = uicontrol('style','push','units','normal','pos',[.92 .79 .07 .06], ...
        'string','RFC','call','democc_ccrfc=democc_rfcdef(democc_y,democc_tp,democc_k,democc_ccrfc); democc_plotmat(democc_demow2,democc_ccrfc,democc_ccmM)');
democc_MM = uicontrol('style','push','units','normal','pos',[.92 .72 .07 .06], ...
        'string','mM','call','democc_ccmM=democc_mmdef(democc_y,democc_tp,democc_k,democc_ccmM); democc_plotmat(democc_demow2,democc_ccrfc,democc_ccmM)');
democc_TP = uicontrol('style','push','units','normal','pos',[.92 .65 .07 .06], ...
        'string','TP','call','democc_tpdef(democc_y)');
democc_TP = uicontrol('style','push','units','normal','pos',[.92 .58 .07 .06], ...
        'string','Redraw','call','democc_k=democc_markmax(democc_y,democc_tp,democc_k,0);');
democc_TP = uicontrol('style','push','units','normal','pos',[.92 .51 .07 .06], ...
        'string','END','call','delete(democc_demow1), delete(democc_demow2), clear democc_*');


