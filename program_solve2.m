% asdf
if !exist('A','var')
  load('program.mat');
  [x, l] = program_init2(A, b, c, v, w, o);
  % current optimal values for primal and dual programs
  op = x'*c * 1.1;
  od = l'*b / 1.1;

  delta = 1/(13*sqrt(m)); % guaranteed to work, via Renegar
  stats = [];
endif

step = 1;
oratio0 = Inf;

while step < 30
  t = time();
  iters = 0;
  [x, iter] = goto_center_primal(A, b, c, op, x);
  iters += iter;
  [l, iter] = goto_center_dual  (A, b, c, od, l);
  iters += iter;
  op2 = c'*x;
  od2 = l'*b;
  op = (1 - delta)*op + delta*op2;
  od = (1 - delta)*od + delta*od2;
  [x2, iter] = goto_center_primal(A, b, c, op, x);
  iters += iter;
  [l2, iter] = goto_center_dual  (A, b, c, od, l);
  iters += iter;
  dx = x - x2;
  dl = l - l2;

  if isinf(oratio0)
    k = 0.75;
  elseif oratio0 > 2
    k = (1 - 1/oratio0);
  else
    k = 0.75;
  endif
  %kappa = 0.3;

  if norm(dx)/norm(x) > eps
    wp = w_to_border_primal(A, b, dx, x2, c, op);
    x = x2 - wp*dx*k;
    s = A*x - b;
    op = c'*x + min([min(x), min(s)]) * 0.5;
    %op2 = c'*(x + dx);
    %op = min(op, op2);
    %c'*[x2,x]
  else
    x = x2;
  endif

  if norm(dl)/norm(l) > 1e-10
    wd = w_to_border_dual  (A, c, dl, l2, b, od);
    l = l2 + wd*dl*k;
    s = c - (l'*A)';
    od = b'*l - min([min(l), min(s)]) * 0.5;
    %od = b'*(l - dl);
  else
    l = l2;
  endif


  dt = time() - t;
  o = [l'*b, x'*c];
  oratio = o(2)/o(1);
  s = [o(1), o(2), dt, iters];
  stats = [stats; s];
  step += 1;

  printf("gap=%f, gapratio=%f\n", oratio - 1, oratio0/oratio - 1);
  if oratio0/oratio - 1 < 1e-6 || oratio - 1 < 1e-2
    break;
  endif

  oratio0 = oratio;
endwhile
sum(stats(:,3:4))
f = fopen(sprintf('res%09i.csv', v), 'w');
fprintf(f, '%i,%g,%i,%g,%g\n', v, sum(stats(:,3)), sum(stats(:,4)), b'*l, c'*x);
fclose(f);

% beter sig som O(v^1.16)
