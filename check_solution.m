% asdf
function [res, xx, oo] = check_solution(A, b, c, x)
  % https://www.matem.unam.mx/~omar/math340/duality.html
  % https://sites.math.washington.edu/~burke/crs/407/notes/section4.pdf
  o = c'*x;
  s = b - A*x;

  res = false;
  xx = x;
  oo = o;

  % construct basis from constraints with smallest slack
  [ss, ssi] = sort(s);
  si = ssi(1:size(A,2));

  if sprank(A(si,:)) < size(A,2) - 1
    printf("bad rank\n");
    return
  endif

  xx = A(si,:)\b(si);
  oo = c'*x;
  % remaining slack
  s2 = b - A*xx;

  if oo <= o-o*eps || min(s2) < -norm(s2)*eps
    o
    oo
    min(s2)
    printf("bad oo or s2\n");
    return
  endif

  % compute dual solution from basis
  y = c'/A(si,:);

  % it suffices to check that y_B >= 0. y_N can be taken as zero
  if min(y) < -norm(y)*eps
    min(y)
    printf("bad y\n");
    return
  endif

  % optionally we could check if the dual solution is close to the primal solution
  % but we don't have enough accuracy to do that
  % do a coarse check anyway
  eps2 = 1e-6;
  ooo = y*b(si);
  if ooo < oo-oo*eps2 || ooo > oo+oo*eps2
    oo
    ooo
    printf("bad ooo\n");
    return
  endif

  res = true;
endfunction
