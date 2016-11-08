function test_suite=test_fitgenpar_mld()
  initTestSuite;
end
function test_fitgenpar_mld_()
    data0 = rndgenpar(0.3,1,0,200,1); 
   x_ML = fzero(@(p)fitgenparml(p, data0), 0); 
   [f,k_ML,s_ML] = fitgenparml(x_ML, data0); % Estimates k_ML and s_ML 
   data1 = floor(data0*10)/10; 
   x = (0:0.1:(max(data1)+0.1))'; 
   N = histc(data1+0.05,x); 
   x_MLD = fzero(@(p)fitgenpar_mld(p,[x, N(:)]) ,0); 
   [f,k_MLD,s_MLD] = fitgenpar_mld(x_MLD,[x N(:)]); % Estimates k_ML and s_ML
end
