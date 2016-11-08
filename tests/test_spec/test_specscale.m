function test_suite=test_specscale()
  initTestSuite;
end
function test_specscale_()
  
   Hm0 = 0.133; Tp = 1.36; 
   Sj = jonswap(linspace(0,125,1025),[Hm0,Tp,3]);  
   ch = spec2char(Sj,{'Hm0','Tm02','Ss'}); 
       % to the corresponding spectrum with Hm0=12 and Ss=ch(3) 
   Ss=ch(3);Tm02=ch(2);Hm0b=12; 
   m0n = (Hm0b/4)^2;  
   m2n = 4*pi^2*m0n*Hm0/(Tm02^2*Hm0b);  
   Sn = specscale(Sj,m0n,m2n,1); 
   ch2 = spec2char(Sn,{'Hm0','Tm02','Ss'}); 
   assert(ch2, [12.0003689455427818, 9.9702039670038527, 0.0773502064191521], 1e-5)
end
