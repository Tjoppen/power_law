% asdf
function w = w_to_border_dual (A, c, dl, l)
  epsilon = 1e-6;
  num = c - (l'*A)';
  den = -(dl'*A)';
  %li = find(abs(den) > epsilon);
  l = num ./ den;
  %l = l(li);
  l = l(find(l > epsilon));
  if length(l) == 0
    error("can't go anywhere");
  endif

  w = min(l);
  if w <= 0
    error("w <= 0");
  endif
endfunction
