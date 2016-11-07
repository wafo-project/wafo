function y = humps(x)
%    Computes a function that has three roots, and some humps.
%
%  Example:
%    x = linspace(0,1);
%    y = humps(x);
%    h = plot(x,y);
%       
%    close();
  if nargin<1 
     x = linspace(0, 1);
  end

  y = 1.0 ./ ((x - 0.3) .^ 2 + 0.01) + 1.0 ./ ((x - 0.9) .^ 2 + 0.04) + 2 .* x - 5.2;
end