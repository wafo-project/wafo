function y = memfun(x,P1,P2,P3,P4)
% MEMFUN private function for MEM  
y = ((simpson(P1,P2.*repmat(exp(P3*transpose(x)),1,P4))));

% comment out this if fmins,fminu, fminsearch fminunc is used
%return

ind = isnan(y);
if any(ind)
  y(ind) = x(ind)*100;
end
%(y)
%y = log(1+y);

%
y = sum((abs(y)));

%