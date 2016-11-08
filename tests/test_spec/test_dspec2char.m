function test_suite=test_dspec2char()
  initTestSuite;
end
function test_dspec2char_()
     S      = demospec('dir'); 
    [ch,txt] = dspec2char(S,1:26);        % fact a vector of integers 
    assert(txt, {'FMdir', 'FPdir','FSpr','FSkew','FMSpr','FLcrst','FS1',... 
       'FS2','FD1','FD2','TpMdir','TpSpr','TpSkew','Wdir','Wdir2','Mdir',... 
       'Pdir', 'Spr', 'Skew', 'MSpr', 'Lcrst', 'S1','S2','D1','D2','TMdir'}) 
 
    ch0 = cell2struct(ch,txt,2);          % Make a structure      
    plot(S.w,ch0.FMdir) 
    assert(dspec2char(S,'wdir'){1}, 0, 1e-10)  % fact a string 
    assert(all([dspec2char(S,{'mdir','pdir'}){:}]<1e-10)) % fact a cellarray of strings 
    assert(all([dspec2char(S,'mdir','pdir'){:}]<1e-10))  % strings 
     
    close all
end
