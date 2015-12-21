% JAPANSEA  coastline map of The Japan Sea
%   
%  CALL:  map = load('japansea.dat');
%  
% Size             :     692 X 2
% Sampling Rate    :    
% Device           :    
% Source           :    http://crusty.er.usgs.gov/coast/getcoast.html
% Format           :    ascii, c1: longitude c2: latitude 
% Description      :
%   JAPANSEA.DAT contains data for plotting a map of The Japan Sea.
%   The data is obtained from USGS coastline extractor.
%   
%  Example: the map is seen by
%
%   plot(map(:,1),map(:,2)), hold on
%   text(131,46,'China')
%   text(132,43,'Vladivostok')
%   text(132,40,'Japan Sea')
%   text(135,35,'Japan')
%   text(139.5,38.3,'Yura')
%   text(139,35.7,'Tokyo'), hold off
%
% % If you have the m_map toolbox (see http://www.ocgy.ubc.ca/~rich/):
%  m_proj('lambert','long',[130 148],'lat',[30 48]);
%  m_line(map(:,1),map(:,2));
%  m_grid('box','fancy','tickdir','out');
%  m_text(131,46,'China');
%  m_text(132,43,'Vladivostok');
%  m_text(132,40,'Japan Sea');
%  m_text(135,35,'Japan');
%  m_text(139.5,38.3,'Yura');
%  m_text(139,35.7,'Tokyo');
% 

% History:
% revised pab Nov 2004
% -Enhanced the example with commands to the m_map toolbox.
% By pab 15.01.2000





