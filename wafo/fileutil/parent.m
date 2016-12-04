function list = parent(h)
%PARENT Lists parents of a Custom object class
%
%  CALL list = parent(h)
%
% 



builtInClasses = {'double','single','logical','char','cell','function_handle','int8',...
      'uint8','int16','uint16','int32','uint32','int64','uint64'};    
 
    fn = fieldnames(h);
    h  = struct(h);
    c  = struct2cell(h(1));
    t  = cellfun(@class,c,'UniformOutput',false);


list = t(~ismember(t,builtInClasses) & strcmp(t,fn));

%list = t(~ismember(t,builtInClasses) ) ;