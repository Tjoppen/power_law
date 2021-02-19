% Goes to the center of the polytope defined by A and b, starting at x0
function [x, iters] = goto_center (A, b, x0)
  x = x0;
  k = 1;
  iters = 0;
  % these are the state used by L-BFGS
  % currently I only use a single iteration's worth of history
  % using more than that seems to be detrimental, probably because this algorithm takes quite large steps
  ss = [];
  ys = [];
  rho = [];

  function y = pr2 (x, P, ss, ys, rho, k)
    if k == 0
      y = P\x;
    else
      z1 = x - ys(:,k)*(rho(k)*(ss(:,k)'*x));
      z2 = pr2(z1, P, ss, ys, rho, k-1);
      z3 = z2 - ss(:,k)*(rho(k)*(ys(:,k)'*z2));
      y = z3 + ss(:,k)*(rho(k)*(ss(:,k)'*x));
    endif
  endfunction

  % preconditioner
  function y = pr (x, P, ss, ys, rho)
    % Used to be the two-loop version of L-BFGS but that didn't work for some reason
    % Do recursive version instead
    y = pr2(x, P, ss, ys, rho, size(ss,2));
  endfunction

  while true
    % compute slacks, sanity check
    s = b - A*x;
    if min(s) <= 0
      [w,iw] = min(s)
      k
      error("min(s) <= 0");
    endif

    s = 1./s;
    g = (s'*A)';
    s2 = s.^2;
    S = diag(s2);

    if k > 1
      % gather stats for L-BFGS
      ys = g - glast;
      rho = [1/(ys'*ss)];
    endif

    % compute diagonal preconditioner
    bsz = 100000;
    p = zeros(1,size(A,2));
    for m = 1:bsz:size(A,2)
      l = min(m+bsz-1, size(A,2));
      p(m:l) = sum(S*A(:,m:l).^2,1);
    endfor
    P = diag(p);
    Pfun = @(x) pr (x, P, ss, ys, rho);

    [h, flags, relres, iter] = pcg(@(x) (((S*(A*x))')*A)', -g, 1e-6, length(b), Pfun, [], zeros(size(x)));
    iters += iter;

    if flags == 4
      error("P is not SPD");
    elseif flags > 1 % accept PCG not reaching the specified tolerance
      flags
      error("solve failed");
    endif

    % we're not always guaranteed to be able to take the full step h
    % try halving h until it works
    rr = 1;
    xn = x + h*rr;
    while min(b - A*xn) <= 0
      rr /= 2;
      xn = x + h*rr;
    endwhile
    x = xn;

    ncur = norm(h)/norm(x);
    gnorm = norm(g);
    printf("iter=%5i gnorm=%g ncur=%g cond≃%g\n", iter, gnorm, ncur, max(p)/min(p));

    if ncur < 1e-6 || gnorm < 1e-10
      return
    endif

    ss = h;
    glast = g;
    k += 1;
  endwhile
endfunction

