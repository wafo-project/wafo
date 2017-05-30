function test_suite=test_parsemfilestr()
  initTestSuite;
end
function test_parsemfilestr_()
           % the dewhite function 
 str =  freadtxt('dewhite');  
 [names,synopsis,subroutines,HH,h1] = parsemfilestr(str,'math','trivia') 
 assert(names, strvcat('cell', 'error','find','iscell','ischar','isempty',... 
       'isspace','max','nargchk','nargin','numel','size'))
end
