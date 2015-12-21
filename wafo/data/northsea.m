% NORTHSEA  coastline map of The Nortsea
%   
%  CALL:  map = load('northsea.dat');
%  
% Size             :     60646 X 2
% Sampling Rate    :    
% Device           :    
% Source           :    http://crusty.er.usgs.gov/coast/getcoast.html
% Format           :    ascii, c1: longitude c2: latitude 
% Description      :
%   NORTHSEA.DAT contains data for plotting a map of The Northsea.
%   The data is obtained from USGS coastline extractor.
%   
%  Example: the map is seen by
%
%   plot(map(:,1),map(:,2)), hold on
%   text( 1,62     ,'Statfjord A'), line([1.8, 1.8], [62  ,61.2 ])
%   plot(1.8,61.2,'x')  
%   text(1,60.5,'Gullfaks C'),      line([1.8, 2.3 ], [60.5 ,61.20 ])
%   plot(2.30,61.20,'x')
%   text(1,59.1,'Frigg'),           line([1.8, 2.0 ], [59.1 ,59.9 ])
%   plot(2.0,59.90,'x'),
%   text(1,57.6,'Sleipner'),        line([1.8, 1.9 ], [57.60 ,58.4 ])
%   plot(1.90,58.40,'x') 
%   text(1,56.9,'Draupner'),        line([1.8, 2.6 ], [56.90 ,57.7 ])
%   plot(2.6,57.7,'x'),
%   text(10.40,60.10,'Oslo'),      plot(10.80,59.85,'h')
%   text(10.00,63.05,'Trondheim'), plot(10.80,63.40,'h')
%   text(4.00,58.80,'Stavanger'),  plot(5.52,58.90,'h')
%   text(3.50,60.30,'Bergen'),     plot(5.20,60.30,'h') , hold off
%
% % If you have the m_map toolbox (see http://www.ocgy.ubc.ca/~rich/):
%  m_proj('lambert','long',[-2 12],'lat',[56 64]);
%  m_line(map(:,1),map(:,2));
%  m_grid('box','fancy','tickdir','out');
%  m_text( 1,62     ,'Statfjord A') ;
%  m_line([1.8, 1.8], [62  ,61.2 ]);
%  m_text(1.7,61.2,'x') ;
%  m_text(1,59.5,'Gullfaks C') ;
%  m_line([1.8, 2.3 ], [59.50 ,61.20 ]);
%  m_text(2.20,61.20,'x') ;
%  m_text(10.10,60.05,'Oslo'); 
%  m_text(8.5,63.5,'Trondheim'); 
%  m_text(4.00,58.80,'Stavanger'); 
%  m_text(3.50,60.30,'Bergen') ;
%  m_text(8,61,'Norway');

%History
% revised pab Nov 2004
% -Enhanced the example with commands to the m_map toolbox
% By pab 06.07.00