function self = label(self)
%LABEL Set automatic x-,y- and z- labels on SPECDATA object
Nf = numel(self);
if Nf>1
  for ix=1:Nf
    self(ix) = label(self(ix));
  end
  return
end


N = length(self.type);

if N==0
   error('Object does not appear to be initialized, it is empty!')
end


if self.norm~=0
  
  switch lower(self.type(N-2:N))
    case 'dir'
      title1 = 'Normalized Directional Spectrum';
      if strcmpi(self.freqtype,'w')
        labels = {'Normalized Frequency' ,'','S(wn,\theta)'};
      else
        labels = {'Normalized Frequency' ,'','S(fn,\theta)'};
      end
      switch lower(self.angletype(1))
        case 'r'
          labels{2} = 'Wave directions [rad]';
        case 'd'
          labels{2} = 'Wave directions [deg]';
      end
    case 'req'
      title1 = 'Normalized Spectral density';
      if strcmpi(self.freqtype,'w')
        labels = {'Normalized Frequency' ,'S(wn)'};
      else
        labels = {'Normalized Frequency' ,'S(fn)'};
      end
    case 'k1d'
      title1 = 'Normalized Spectral density';
      labels = {'Normalized Wave number','S(kn)'};
    case 'k2d'
      title1 = 'Normalized Wave Number Spectrum';
      xlbl_txt = 'Normalized Wave number';
      labels = {xlbl_txt , xlbl_txt,'S(kn1,kn2)'};
    otherwise
       error('Object does not appear to be initialized, it is empty!')
     
  end
   
else
  
 switch lower(self.type(N-2:N))
    case 'dir'
      title1 = 'Directional Spectrum';
      if strcmpi(self.freqtype,'w')
        labels = {'Frequency [rad/s]' ,'','S(w,\theta) [m^2 s / rad^2]'};
      else
        labels = {'Frequency [Hz]' ,'','S(f,\theta) [m^2 s / rad]'};
      end
      switch self.angletype(1)
        case 'r'
          labels{2} = 'Wave directions [rad]';
        case 'd'
          labels{2} = 'Wave directions [deg]';
      end
    case 'req'
      title1 = 'Spectral density';
       if strcmpi(self.freqtype,'w')
        labels = {'Frequency [rad/s]' ,'S(w) [m^2 s / rad]'};
      else
        labels = {'Frequency [Hz]' ,'S(f) [m^2 s]'};
      end
      
    case 'k1d'
      title1 = 'Spectral density';
      labels = {'Wave number [rad/m]','S(k) [m^3/ rad]'};
    case 'k2d'
      title1 = 'Wave Number Spectrum';
      xlbl_txt = 'Wave number [rad/m]';
      labels = {xlbl_txt , xlbl_txt,'S(k1,k2) [m^4/ rad^2]'};
   otherwise
      error('Object does not appear to be initialized, it is empty!')
 end
end
self.wdata = set(self.wdata,'labels',labels,'title',title1);   
