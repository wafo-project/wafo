function args = _mesh_args(f)
  if length(f.args{1})==numel(f.args{1})
     n = length(f.args);
     args = cell(1,n);
     [args{1:n}] = meshgrid(f.args{:});
   else
     args = f.args;
   end
end