% like goto_center() except c and eye(n) in A are implicit
function [x, iters] = goto_center_primal(A, b, c, op, x0)
  x = x0;
  iters = 0;

  function y = fun (z, x2, d, S)
    if any(isnan(z))
      error('NaN')
    endif
    if any(isinf(z))
      error('Inf')
    endif
    %printf("%.1f\t", log10(norm(z)));
    printf(".");
    fflush(stdout);
    y1 = (((S*(A*z))')*A)';
    y2 = x2.*z;
    y3 = c*((c'*z)*d^-2);
    y = y1 + y2 + y3;
  endfunction

  while true
    s = A*x - b;
    d = op - c'*x;
    if min(s) <= 0 || d <= 0
      min(s)
      d
      error('bad s or d');
    endif

    g = ((1./s)'*A)' + 1./x - c/d;
    S = diag(s.^-2);
    x2 = x.^-2;

    % compute diagonal preconditioner
    bsz = 100000;
    p = zeros(1,size(A,2));
    for m = 1:bsz:size(A,2)
      printf(",");
      fflush(stdout);
      l = min(m+bsz-1, size(A,2));
      p(m:l) = sum(S*A(:,m:l).^2,1);
    endfor
    p += x'.^-2;
    p += c'.^2/d^2;
    P = diag(p);

    [h, flags, relres, iter] = pcg(@(z) fun(z, x2, d, S), -g, 1e-6, 100, P, [], -(P\g));
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
    xn = x - h*rr;
    while min(A*xn - b) <= 0 || op - c'*xn <= 0 || min(xn) <= 0
      printf("r");
      fflush(stdout);
      rr /= 2;
      xn = x - h*rr;
    endwhile
    x = xn;

    nh = norm(h)/norm(x);
    printf(" prim norm(h)=%g iter=%i\n", nh, iter);
    if nh < 1e-10
      return
    endif
  endwhile
endfunction
