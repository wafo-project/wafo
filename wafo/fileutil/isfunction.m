function val = isfunction(file)
%ISFUNCTION Return true if m-file is a function.


  if iscell(file)
    val = zeros(size(file));
    for ix = 1:length(file)
      try
        _nargin(file{ix});
        val(ix) = 1;
      end
    end
  else
    val = 0;
    try
      _nargin(file);
      val = 1;
    end
  end
end

function nout = _nargin(file)
  try
    nout = nargin(file)
   catch
     [d, n, e] = fileparts(file);
     nout = nargin(n);
   end
end