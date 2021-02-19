% like goto_center_primal() except for the dual program
function [l, iters] = goto_center_dual(A, b, c, od, l0)
  l = l0;
  iters = 0;

  function y = fun (z, l2, d, S)
    if any(isnan(z))
      error('NaN')
    endif
    if any(isinf(z))
      error('Inf')
    endif
    %printf("%.1f\t", log10(norm(z)));
    printf(".");
    fflush(stdout);
    y1 = A*(S*(z'*A)');
    y2 = l2.*z;
    y3 = b*((b'*z)*d^-2);
    y = y1 + y2 + y3;
  endfunction

  nh0 = 0;

  while true
    s = c - (l'*A)';
    d = b'*l - od;
    if min(s) <= 0 || d <= 0
      min(s)
      d
      error('bad s or d');
    endif

    g = A*(1./s) + 1./l - b/d;
    S = diag(s.^-2);
    l2 = l.^-2;

    % compute diagonal preconditioner
    bsz = 100000;
    p = zeros(1,size(A,1));
    for m = 1:bsz:size(A,1)
      printf(",");
      fflush(stdout);
      o = min(m+bsz-1, size(A,1));
      p(m:o) = sum(A(m:o,:).^2*S,2);
    endfor
    p += l'.^-2;
    p += b'.^2/d^2;
    P = diag(p);

    [h, flags, relres, iter] = pcg(@(z) fun(z, l2, d, S), -g, 1e-6, length(c), P, [], -(P\g));
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
    ln = l - h*rr;
    while min(c' - ln'*A) <= 0 || b'*ln - od <= 0 || min(ln) <= 0
      printf("r");
      fflush(stdout);
      rr /= 2;
      ln = l - h*rr;
    endwhile
    l = ln;

    nh = norm(h)/norm(l);
    printf(" dual norm(h)=%g iter=%i\n", nh, iter);
    if nh < 1e-10
      return
    endif
  endwhile
endfunction
