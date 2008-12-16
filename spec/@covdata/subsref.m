function b = subsref(s,index)
% SUBSREF Define field name indexing for SPECDATA objects
%
% 

N = length(index);
b = s;
for ix = 1:N
  subs = index(ix).subs;
  switch index(ix).type
    case '()'
      b = b(subs{:});
    case '.'
      fn   = fieldnames(b);     
      
      if any(strcmp(subs,fn))
        b = b.(subs);
      else
        switch lower(subs)
          case {'getcovariance', 'r','data'}
            b = get(b.wdata,'data');
          case {'x','y','t'}
            x = get(b.wdata,'args');
            isND = length(b.lagtype)>1;
            dim = find(b.lagtype==subs); 
            if isempty(dim)
              error('???? Reference to non-existent field or method %s',subs)
            end

            if isND ~= iscell(x)
              error('args of WDATA object must be a cellarray for 2D or 3D data and vector for 1D data!')
            end
            if isND 
              b = x{dim};
            else
              b = x;
            end
            
          otherwise
            %case {'name','type','data','args','argMin','argMax', 'sampleRate', 'labels',    'title',...
            %    'contourLevels', 'percentLevels',  'workspace',    'info',   'note', 'date'}
            
            try
              %b = get(s.wdata,subs);
              b = subsref(b.wdata,index(ix:end));
              return
            catch
              error('???? Reference to non-existent field or method %s',subs)
            end
        end
      end
  end
  
end

