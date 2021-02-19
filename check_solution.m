% Finds the closest extreme point to x and checks whether that is a valid solution
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
  k = sprank(A(si,:))

  if k != size(A,2)
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
    printf("objective worse, or basis is not feasible\n");
    return
  endif

  % compute dual solution from basis
  y = c'/A(si,:);

  % it suffices to check that y_B >= 0, I think. y_N can be taken as zero
  if min(y) < -norm(y)*eps
    min(y)
    printf("dual not feasible\n");
    return
  endif

  % optionally we could check if the dual solution is the same as the primal solution
  eps2 = 1e-2; % being within 1% of optimal is good enough for us
  ooo = y*b(si);
  relo = (ooo - oo) / oo
  if relo > eps2
    printf("not quite there yet\n");
    return
  endif

  res = true;
endfunction
