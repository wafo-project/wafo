function test_suite=test_locate()
  initTestSuite;
end
function test_locate_()
  
   loc0 = locate('*plot.m');    % Find all files ending with 'plot.m'. 
   loc1 = locate('im*');        % Find all files starting with 'im'. 
   names = {}; 
   for i=1:length(loc0), 
    [folder, names{i}] = fileparts(loc0{i}); 
   end 
   i0 = [strmatch('pdfplot', names)(1), strmatch('trplot', names)(1), ... 
         strmatch('nplot', names)(1)];  
   assert(names(i0), {'pdfplot','trplot','nplot'});
end
