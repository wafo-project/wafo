function [pvalue, test, v] = chi2gof(fs,ft)
% CHI2GOF CHI Squared Goodness Of Fit test.
%
% CALL:  P = chi2gof(fs, ft)
%
%   P      = P - value  
%   fs     = fitted PDF structure/matrix
%            evaluated at observed points.
%   ft     = theoretical PDF structure/matrix
%            evaluated on a equidistant grid 
%
%  Large P value -> good fit
%  Small P value -> lesser fit
%
%  Example: % Check how well rayleigh data can be described by N(0,1)
%   xs = rndray(1,500,1);
%   x  = linspace(-7,7,201); 
%   p  = chi2gof(pdfnorm(xs),pdfnorm(x));
%           
% See also  cdfchi2, qlevels


%Tested on: matlab 5.3
% History:
% revised pab 11.11.2000
% - made it independent of stats toolbox only dependent on wstats
% by pab 21.09.99  
  


% CHI^2 goodness of fit (GOF)   test

 if isstruct(fs) % structure
     r2=fs.f(:);
   else
    r2 = fs(:);
  end
ntresh=length(r2);
if true
  %ntresh>120 % only valid for samples larger than 120
  if ntresh<30
    warning('WAFO:CHI2GOF','Too few data -> interpret result with caution!')
  end

   if isstruct(ft)
     r=ft.f;
   else
    r = ft;
  end
  k=max([ceil(sqrt(ntresh)),8]);%divide the data into k subsets (greater than 8)
  pk=100/k;                     % with equal probabilty
    
  grpdiv=qlevels(r,(100-pk):-pk:pk); % find dividing levels 
  %grpdiv=pdfnorm(norminv( (pk:pk:(100-pk))/200))' 
  
  np=ntresh/k; %the expected number of points in each group (must be greater than 5)
  
  Nobs = histc(r2,[-inf;grpdiv(:);inf]);
  Nobs(k+1:end) = [];
  if np<5 || any(abs(Nobs-np)>np)
    %GTEST
    test = 2*sum(Nobs.*log(Nobs/np+realmin))+2*sum(Nobs==0)*log(realmax);
  else
    %Pearson's chi-squared test statistic
    
    test = sum((Nobs-np).^2/np);
    %     test=(sum(grpdiv(1 )>=r2)-np)^2/np +(sum(grpdiv(k-1)<r2)-np)^2/np;
    %     for ix=2:k-1,
    %       test=test+(sum((grpdiv(ix-1 )<r2).*(grpdiv(ix )>=r2))-np)^2/np;
    %     end
    
  end
  %test %test statistic
  pvalue=cdfchi2(test,k-1,'lowertail',false); % pvalue
  v=k-1;
else  
  pvalue=[];
  disp('to few data')
  %ntresh
end
return


