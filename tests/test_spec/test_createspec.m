function test_suite=test_createspec()
  initTestSuite;
end
function test_createspec_()
     S=createspec('dir'); 
    S.date = ''; 
    names = {'S', 'w', 'theta', 'tr', 'h', 'type', 'phi','norm','note', 'date'}; 
    vals = {[],[],[],[],inf,'dir', 0, 0, [], ''}; 
    assert(fieldnames(S), names'); 
    assert(struct2cell(S), vals');
end
