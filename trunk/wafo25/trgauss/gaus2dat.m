function x = gaus2dat(xx,g)
% GAUS2DAT Transforms  xx  using the inverse of transformation  g.
%
%  CALL: x = gaus2dat(xx,g);
%
%        xx = input data with time in first column and values in second.
%        g  = transform function, xx=g(x).
%        x  = invers transformed data, x=G(xx), G=g^(-1).
%
% See also  dat2gaus  and  tranproc.

x=[xx(:,1), tranproc(xx(:,2:size(xx,2)),fliplr(g))];
