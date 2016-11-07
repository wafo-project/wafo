function out = split_linespec(linespec, colorname)
% Split linespec into linestyle, marker and color
%
%  Example
%  assert(split_linespec('b'), struct('color', 'b'));
%  assert(split_linespec('b--'), struct('linestyle','--', 'color', 'b'));
%  assert(split_linespec('b.'), struct('marker', '.', 'color', 'b'));
%  assert(split_linespec('b:^'), struct('linestyle',':', 'marker', '^', 'color', 'b'));
%  assert(split_linespec('blablar'), struct('color', 'r'));
%  assert(split_linespec('b', 'linecolor'), struct('linecolor', 'b'));

  if nargin<2 || isempty(colorname)
    colorname = 'color';
  end
  line_styles = {'--', '-.','-', ':'};
  markers = {'\+','o','\*','\.','x','s','square','d','diamond','\^','v','>','<',...
      'p', 'pentagram', 'h', 'hexagram'};
  colors = {'r','g','b', 'c', 'm','y','k','w'};
  out = struct();
  [linestyle, noMatch] = get_style(linespec, line_styles);
  [marker, noMatch] = get_style(noMatch, markers);
  color = get_style(noMatch, colors);
  if length(linestyle), out.linestyle = linestyle; end
  if length(marker), out.marker = marker; end
  if length(color), out.(colorname) = color; end
end

function [out, remainder] = get_style(linespec, valid_styles)
  remainder = linespec;
  out = '';
  if ~any(linespec)
    return
  end
  for style = valid_styles,
    [match, noMatch] = regexp(linespec, sprintf('(%s){1}',style{1}), 'match', 'split');
    if length(match)==1,
      out = match{1};  % return only first occurrence
      remainder = strcat(noMatch{:});
      return
    end
  end
end
