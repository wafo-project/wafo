function test_suite=test_specdata()
  initTestSuite;
end
function test_specdata_()
  Sold = createspec; 
 fn   = fieldnames(Sold); 
 H = mkjonswap; 
 w = linspace(0,4,256)'; 
 S = specdata(H(w),w);  % Make spectrum object from numerical values 
  
 Sold = get(S,fn{:});    % Convert to old type of spectrum object.
end
