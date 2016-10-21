function s = subsasgn(s,index,b)
% SUBSASGN Define field name assigment for SPECDATA objects
%
% 

switch index.type
case '()'
   s(index.subs) = b;
  case '.'
    fn   = fieldnames(s);
    subs = index.subs;
    if ismember(subs,fn)
      if strcmpi(subs,'wdata') && ~isa(b,'wdata')
        error('Not allowed to assign a non-WAFO data object to wdata.')
      end
      s.(index.subs) = b;
    else
       switch lower(subs)
         case {'setspectrum', 's','data'}
         case {'f','w','theta','k','k2'}
           x = get(s.wdata,'args');
           is2D = ismember(s.type,{'k2d','dir','encdir'});
           if is2D ~= iscell(x)
             error('args of WDATA object must be a cellarray for 2D data and vector for 1D data!')
           end
           argMin = get(s.wdata,'argMin');
           argMax = get(s.wdata,'argMax');

           if is2D
             if strcmpi(s.type(end-2:end),'dir')
               dim = 2*strcmpi(subs,'theta') + strcmpi(subs,s.freqtype);
             else
               dim = 2*strcmpi(subs,'k') + strcmpi(subs,'k2');
             end
             if dim ==0
               error('???? Reference to non-existent field or method %s',subs)
             end
             x{dim} =  b;
             argMin(dim) = min(b);
             argMax(dim) = max(b);
           else
             if s.freqtype~=subs
               error('???? Reference to non-existent field or method %s',subs)
             end
             x = b;
             argMin = min(b);
             argMax = max(b);
           end
            s.wdata = set(s.wdata,'args',x,'argMin',argMin,'argMax',argMax);
         otherwise
           %case {'name','type','data','args','argMin','argMax', 'sampleRate', 'labels',    'title',...
           %    'contourLevels', 'percentLevels',  'workspace',    'info',   'note', 'date'}
           
           try
             subsasgn(s.wdata,index,b);
           catch
             error('???? Reference to non-existent field or method %s',subs)
           end
       end
    end
    
    %error('???? Reference to non-existent field or method %s',subs)
end
