function recdemo
% RECDEMO  Show Rec figures with point and click interface
%
% CALL:    recdemo
%

% By pab 28.01.2000

header = 'Reconstruction Demo';
labels =strvcat('Change Settings' , 'Run All Scripts' , ...
          'Run One Script' , 'Clean Up');
callbacks = str2mat('recfig(-1)', 'recfig(-2)', ...
             'recfig(-3)', 'recleanup');
global RECFIGNUM 
RECFIGNUM = 0;
clc;
%figure(gcf+1)
recintro
choices('Rec',header,labels,callbacks);
    
