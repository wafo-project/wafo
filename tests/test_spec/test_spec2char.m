function test_suite=test_spec2char()
  initTestSuite;
end
function test_spec2char_()
    S      = demospec; 
   [ch R,txt] = spec2char(S,[1 2 3]);    % fact a vector of integers 
   ch0 = cell2struct(num2cell(ch),txt,2);          % Make a structure  
   assert(ch, [6.80042799184815, 9.66549688067646, 9.41922471783712], 1e-10) 
   assert(txt, {'Hm0', 'Tm01', 'Tm02'})  
   assert(spec2char(S,'Ss'), 0.0491112394396479, 1e-10)     % fact a string 
   assert(spec2char(S,{'Tp','Tp1'}), ... 
                      [11.0231321178589, 10.8684288697414], 1e-10)   % fact a cellarray of strings
end
