function test_suite=test_sensortype()
  initTestSuite;
end
function test_sensortype_()
  validNames = {'n','n_t','n_tt','n_x','n_y','n_xx','n_yy','n_xy','p','u',... 
               'v','w','u_t', 'v_t','w_t','x_p','y_p','z_p','nan'}; 
 for i=1:19, 
   assert(sensortype(i), validNames{i}); 
 end
end
