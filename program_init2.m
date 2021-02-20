% asdf
function [x, l] = program_init2(A, b, c, v, w, o)
  [m, n] = size(A);
  k = 1.1;

  x = bicgstab(A(1:v,1:v), b(1:v));
  x = k*x./min((A(1:(v+w),1:v)*x)./b(1:(v+w)));

  x((v+1):n,:) = -A((v+w+1):m,1:v)*x(1:v,:) * k;

  if min(x) <= 0 || min(A*x - b) <= 0
    min(x)
    min(A*x - b)
    error('bad x');
  endif

  %l = [1./b(1:(v+w)); 0.9*ones(o,1)];
  l = [1./b(1:(v+w)); 0.9e-6*ones(o,1)];
  %l = ((c-0.1)'/A)';
  %l = (ones(1,v)/A(1:v,1:v))';
  %l = [l; 0.9e-6*ones(o,1)];
  %l = [l./b(1:v);0.9e-6*ones(w+o,1)];

  if min(l) <= 0 || min(c'-l'*A) <= 0
    % HACKHACK
    l((v+w+1):end) *= 1.1;
    if min(l) <= 0 || min(c'-l'*A) <= 0
      error('bad l');
    endif
  endif

  l'*b
  c'*x
endfunction
