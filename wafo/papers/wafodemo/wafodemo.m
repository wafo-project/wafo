function wafodemo
% WAFODEMO  Show WAFO figures with point and click interface
%
% CALL:    wafodemo
%

% tested on: Matlab 5.2
% History
% By pab 28.01.2000


% TODO % update calls to spec2XXpdf programmes due to change in call

header = 'WAFO Demo';
labels =strvcat('Change Settings' , 'Run All Scripts' , ...
          'Run One Script' , 'Clean Up');
callbacks = str2mat('wafofig(-1)', 'wafofig(-2)', ...
             'wafofig(-3)', 'wafocleanup');
global WAFOFIGNUM 
WAFOFIGNUM = 0;
clc; 
wafointro
choices('WAFO',header,labels,callbacks);
    
