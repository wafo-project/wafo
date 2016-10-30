function Gf = peakEnhancementFactor(wn,options)
  %PEAKENHANCEMENTFACTOR 
  %
  % 

  gam = options.gamma;
  sb  = options.sigmaB;
  sa  = options.sigmaA;


  k    = (wn>1);
  sab      = zeros(size(wn));
  sab(k)   = sb;
  sab(~k)  = sa;
  wn(wn<0) = 0;

  wnm12 = 0.5*((wn-1)./sab).^2;
  Gf    = gam.^(exp(-wnm12));

end % peak enhancement factor
