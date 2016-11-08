function test_suite=test_psi()
  initTestSuite;
end
function test_psi_()
   assert(psi([1,2,4,8]), [-0.577215664901553, 0.422784335098447, ... 
									1.256117668431780, 2.015641477955590], 1e-12)
end
