function xx = dat2gaus(x,g)
% DAT2GAUS Transforms  x  using the transformation  g.
%
%  CALL: xx = dat2gaus(x,g);
%
%        x  = input data with time in first column and values in second.
%        g  = transform function, xx=g(x).
%        xx = transformed data, xx=g(x).
%
%  See also  gaus2dat  and  tranproc.

xx=[x(:,1), tranproc(x(:,2:size(x,2)),g)];
