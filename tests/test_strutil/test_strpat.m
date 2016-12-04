function test_suite=test_strpat()
  initTestSuite;
end
function test_strpat_()
  
	assert(strrep('aa_a_aaa_aa','a','XXX'), 'XXXXXX_XXX_XXXXXXXXX_XXXXXX') 
	assert(strpat('aa_a_aaa_aa','a','XXX'), 'aa_XXX_aaa_aa') 
	assert(strpat('aa_a_aaa_aa','aa','X'),'X_a_aaa_X') 
	assert(strpat('aa_a_aaa_aa','a','XXX'),'aa_XXX_aaa_aa') 
  
 %  The following does not work on octave. 
	% assert(strpat([1:3,pi,5:7],pi,nan), [1 2 3 NaN 5 6 7]) 
	% assert(strpat(pi*(1:6),pi*(4:5),[]), [3.1416 6.2832 9.4248 18.85]) 
	% assert(strpat(pi*(1:6),pi*(4:5),[nan inf -inf nan]), ... 
 %              [3.1416 6.2832 9.4248 NaN Inf -Inf NaN 18.85]) 
	% assert(strpat(pi*(1:6),pi*(4:5),1:3), [3.1416 6.2832 9.4248 1 2 3 18.85]) 
	% assert(strpat([-10,nan,1,-10,nan,nan,1,nan,1,-10],[nan,1],inf), ... 
 %              [-10 Inf -10 NaN Inf Inf -10])
end
