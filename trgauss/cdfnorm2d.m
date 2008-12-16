function bvn = cdfnorm2d( b1, b2, r ) 
%CDFNORM2D Bivariate Normal cumulative distribution function
%
%  CALL:  val = cdfnorm2d( b1, b2, r )
%
%    bvn   = distribution function evaluated at b1, b2.  
%    b1,b2 = upper integration limits
%    r     = correlation coefficient  (-1 <= r <= 1).
%
%  CDFNORM2D computes bivariate normal probabilities, i.e.,
%  the probability Prob(X1 <= B1 and X2 <= B2) with an absolute error
%  less than 1e-15.
%     
%  This function is based on the method described by
%  drezner, z and g.o. wesolowsky, (1989),
%  on the computation of the bivariate normal integral,
%  journal of statist. comput. simul. 35, pp. 101-107,
%  with major modifications for double precision, and for |r| close to 1.
%
% Example 
%  x = linspace(-5,5,20);  
%  [B1,B2] = meshgrid(x);
%  r  = 0.3;
%  F = cdfnorm2d(B1,B2,r);  
%  surf(x,x,F)
%
% See also  cdfnorm
  
% Reference
%  Drezner, z and g.o. Wesolowsky, (1989),
%  "on the computation of the bivariate normal integral",
%  journal of statist. comput. simul. 35, pp. 101-107,
  
%History  
% Revised pab 29 march 2007
% -fixed a bug for r = 0.75 thanks to Richard Chiburis
% revised pab 19.01.2001
%  -vectorized the code to handle multiple integration limits  
%  -renamed from bvnd to bvnormcdf  
%  -replaced call to erf with erfc  
%   
% by  alan genz
%     department of mathematics
%     washington state university
%     pullman, wa 99164-3113
%     email : alangenz@wsu.edu
error(nargchk(3,3,nargin))
[csize,h,k,r] = comnsize(-b1,-b2,r);
if any(isnan(csize))
  error ('b1, b2 and r must be of common size or scalar');
end

bvn  = zeros(size(h));

%zero = 0.d0;
%one  = 1.d0;
two  = 2.d0;
twopi = 6.283185307179586d0;
%     gauss legendre points and weights, n = 6
w6 = [ 0.1713244923791705D+00, 0.3607615730481384D+00, 0.4679139345726904D+00];
x6 = [-0.9324695142031522D+00,-0.6612093864662647D+00,-0.2386191860831970D+00];
%     gauss legendre points and weights, n = 12
w12 = [ 0.4717533638651177d-01, 0.1069393259953183d+00, 0.1600783285433464d+00,...
       0.2031674267230659d+00, 0.2334925365383547d+00, 0.2491470458134029d+00];
x12=[ -0.9815606342467191d+00,-0.9041172563704750d+00,-0.7699026741943050d+00, ...
    -0.5873179542866171d+00,-0.3678314989981802d+00,-0.1252334085114692d+00];
%     gauss legendre points and weights, n = 20
w20 = [ 0.1761400713915212d-01, 0.4060142980038694d-01, ...
      0.6267204833410906d-01,  0.8327674157670475d-01, ...
      0.1019301198172404d+00, 0.1181945319615184d+00,...
      0.1316886384491766d+00, 0.1420961093183821d+00,...
      0.1491729864726037d+00, 0.1527533871307259d+00];
x20=[ -0.9931285991850949d+00, -0.9639719272779138d+00, ...
      -0.9122344282513259d+00, -0.8391169718222188d+00, ...
       -0.7463319064601508d+00,-0.6360536807265150d+00,...
       -0.5108670019508271d+00,-0.3737060887154196d+00, ...
   -0.2277858511416451d+00, -0.7652652113349733d-01];

hk  = h.*k;

k0 = find(abs(r) < 0.925d0 );
if (k0) 
  hs = ( h(k0).^2 + k(k0).^2 )/two;
  asr = asin(r(k0));
  k1 = find(r(k0)>=0.75);
  if any(k1)
    k01 = k0(k1);
    for i = 1:10
      for is = -1:2:1, 
        sn = sin( asr(k1).*(is.*x20(i) +1 )/2 );
        bvn(k01) = bvn(k01) + w20(i)*exp( ( sn.*hk(k01) - hs(k1) )./(1 - sn.*sn ) );
      end 
    end 
  end
   k1 = find(0.3 <= r(k0) & r(k0)<0.75);
  if any(k1)
    k01 = k0(k1);
    for i = 1:6
      for is = -1:2:1, 
        sn = sin( asr(k1).*(is.*x12(i) +1 )/2 );
        bvn(k01) = bvn(k01) + w12(i)*exp( ( sn.*hk(k01) - hs(k1) )./(1 - sn.*sn ) );
      end 
    end 
  end
   k1 = find( r(k0)<0.3);
  if any(k1)
    k01 = k0(k1);
    for i = 1:3
      for is = -1:2:1, 
        sn = sin( asr(k1).*(is.*x6(i) +1 )/2 );
        bvn(k01) = bvn(k01) + w6(i)*exp( ( sn.*hk(k01) - hs(k1) )./(1 - sn.*sn ) );
      end 
    end 
  end
  bvn(k0) = bvn(k0).*asr/( two*twopi ) + fi(-h(k0)).*fi(-k(k0));
end

k1 = find(0.925<=abs(r) & abs(r)<=1 );
if any(k1)
  k2 = find(r(k1) < 0);
  if any(k2 ) 
    k12 = k1(k2);
    k(k12)  = -k(k12);
    hk(k12) = -hk(k12);
  end
  k3 = find( abs(r(k1)) < 1);
  if (k3)
    k13 = k1(k3);
    as = (1 - r(k13) ).*(1 + r(k13) );
    a  = sqrt(as);
    b  = abs( h(k13) - k(k13) );
    bs = b.*b;
    c  = ( 4.d0 - hk(k13) )/8.d0;
    d  = ( 12.d0 - hk(k13) )/16.d0;
    asr = -( bs./as + hk(k13) )/2.d0;
    k4 = find(asr > -100.d0);
    if (k4) 
      bvn(k13(k4)) = a(k4).*exp(asr(k4)).*(1 - c(k4).* ....
	  ( bs(k4) - as(k4)).*(1 - d(k4).*bs(k4)/5 )/3 ...
	  + c(k4).*d(k4).*as(k4).^2/5 );
    end
    k5 = find(hk(k13) < 100.d0);
    if ( k5 )
      %               b = sqrt(bs);
      k135 = k13(k5);
      bvn(k135) = bvn(k135) - exp( -hk(k135)/2 ).*sqrt(twopi).*fi(-b(k5)./a(k5)).*b(k5)...
	  .*(1- c(k5).*bs(k5).*(1 - d(k5).*bs(k5)/5 )/3 );
    end
    a = a/two;
    for i = 1:10
      for is = -1:2:1,
	xs  = ( a.*( is*x20(i) + 1 ) ).^2;
	rs  = sqrt(1 - xs );
	asr = -( bs./xs + hk(k13) )/2;
	k6 = find( asr > -100.d0 ) ;
	if any(k6) 
	  k136 = k13(k6);
	  bvn(k136) = bvn(k136) + a(k6).*w20(i).*exp( asr(k6) )...
	      .*( exp(-hk(k136).*(1-rs(k6))./(2*(1+rs(k6))))./rs(k6)...
	      - (1 + c(k6).*xs(k6).*(1 + d(k6).*xs(k6) ) ) );
	end 
      end 
    end 
    bvn(k3) = -bvn(k3)/twopi;
    %bvn = -bvn/twopi;
  end
  k7 = find(r(k1)>0);
  if any(k7 ),
    k17 = k1(k7); 
    bvn(k17) = bvn(k17) + fi( -max( h(k17), k(k17)) );
  end
  k8 = find(r(k1)<0);
  if any(k8), 
    k18 = k1(k8);
    bvn(k18) = -bvn(k18) + max(0, fi(-h(k18)) - fi(-k(k18)) );
  end
end

k9 = find(abs(r)>1);
if any(k9)
  bvn(k9) = NaN;
end
%val = bvn;
return

function F = fi(x)
F = 0.5.*(erfc((-x)./sqrt(2)));
  
