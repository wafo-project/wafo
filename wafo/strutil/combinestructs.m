function s1 = combinestructs(s1,s2,mode)
%COMBINESTRUCTS Combines 2 structures
%
%  CALL s1 = combinestructs(s1,s2,mode)
%
%  s1, s2 = structs to be combined.
%  mode   = Defines the combination mode:           
%           'replace' : replace common fields in s1 with those in s2 (default) 
%           'merge'   : merge all fields of s1 and s2. Common fields take
%                       the values from s2.
%           'add'     : add non-common fields of s2 into s1.
%            
% 
% COMBINESTRUCTS combines the structures S1 and S2 according to selected
% MODE. If MODE is REPLACE the common fields in S2 overwrite the
% corresponding old parameters in S1. If MODE is MERGE the fields in S2 are
% copied to S1. If MODE is ADD the fields in S2 that are not in S1 are
% copied to S1.
%
% Examples
%  s1 = struct('string','hello','lengths',[5 3],'test','tull'); 
%  s2 = struct('string','apple','lengths',[15 30],'color','red'); 
%  s3 = struct('string','apple','lengths',[15 30],'test','tull');
%  s4 = struct('string','apple','lengths',[15 30],'test','tull', 'color', 'red');
%  s5 = struct('string','hello','lengths',[5 3],'test','tull', 'color', 'red');
%  assert(combinestructs(s1,s2,'replace'), s3)
%  assert(combinestructs(s1,s2,'merge'), s4)
%  assert(combinestructs(s1,s2,'add'), s5)
%
% See also parseoptions

if nargin<3 || isempty(mode)
  mode = 'replace';
end


if nargin>=2 && ~isempty(s2)                      % [] is a valid options argument
  N1 = numel(s1);
  N2 = numel(s2);
 
  
  
    s2Names  = fieldnames(s2).';
    s1Names  = fieldnames(s1).';
  
    switch mode
      case {'a','add'} % add
        fnames = setdiff(s2Names,s1Names);
      case {'m','merge'}
        fnames = s2Names;
      otherwise
        % replace common fields
        fnames = intersect(s2Names,s1Names);
        
        
% Old call      
%        [fnames, i2, i1] = intersect(s2Names,s1Names);
%         if isempty(i2),
%           return;
%         end
%         c2 = struct2cell(s2);
%         c1 = struct2cell(s2);
%         c1(i1) = c2(i2);
%         s1 = cell2struct(c1,s1Names,1);

    end
    flag  =  2*(N1>1) + N2>1;
    
    
    N = max(N1,N2);
    
    switch flag
      case 0
        for fn = fnames
          s1.(fn{1}) = s2.(fn{1});
        end
      case 1
        for fn = s1Names
          [s1(1:N).(fn{1})] = deal(s1(1).(fn{1}));
        end
        for fn = fnames
          [s1(1:N).(fn{1})] = deal(s2.(fn{1}));
        end
      case {2,3}
        for fn = fnames
          [s1.(fn{1})] = deal(s2.(fn{1}));
        end
    end
end


return


