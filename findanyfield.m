function [data,cnt] = findanyfield(S,fld)
[data,cnt] = recfun(S,fld,0);
end
function [data,out] = recfun(S,fld,inp)
  data = [];
  out  = -1;
  if isstruct(S)
    if isfield(S,fld)
      data = S.(fld);
      out  = inp+1;
    else
      names = fieldnames(S);
      for k = 1:numel(names)
        [data,out] = recfun(S.(names{k}),fld,inp+1);
        if out>0
          break
        end
      end
    end
  end
end