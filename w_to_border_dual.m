% asdf
function w = w_to_border_dual (A, c, dl, ll, b, od)
  epsilon = 1e-6;
  num = c - (ll'*A)';
  den = -(dl'*A)';
  %li = find(abs(den) > epsilon);
  num = [num; od - b'*ll];
  den = [den; b'*dl];
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

  l = -ll ./ dl;
  l = find(l > epsilon);
  w = min([w; l]);
endfunction
