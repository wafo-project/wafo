function test_suite=test_str2lines()
  initTestSuite;
end
function test_str2lines_()
    s = ['Hello world' char(10) 'Hello space']; 
   lines = str2lines(s); 
   assert(lines, {'Hello world'; 'Hello space'}) 
 
   %To convert from a list of lines back to a string using DOS (CR+LF) line 
   %endings, use 
 
   c = cell(2, length(lines)); 
   c(1,:) = lines; 
   c(2,:) = {sprintf('\n')}; 
   assert(s, [c{1:end-1}]) 
 
   % or 
 
   c = [lines' ; repmat({sprintf('\n')}, 1, length(lines))]; 
   assert(s, [c{1:end-1}])
end
