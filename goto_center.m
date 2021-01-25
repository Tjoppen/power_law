% asdf
%function [x, iters] = goto_center (A, At, A2, b, x0)
function [x, iters] = goto_center (A, b, x0)
  x = x0;
  k = 1;
  iters = 0;
  while true
    s = b - A*x;
    if min(s) <= 0
      [w,iw] = min(s)
      k
      error("min(s) <= 0");
    endif
    s = 1./s;
    g = (s'*A)';
    gnorm = norm(g);
    s2 = s.^2;
    S = diag(s2);
    %SA = S*A;
    %B = At*SA;
    %h = B\-g;
    bsz = 100000;
    p = zeros(1,size(A,2));
    for k = 1:bsz:size(A,2)
      l = min(k+bsz-1, size(A,2));
      p(k:l) = sum(S*A(:,k:l).^2,1);
    endfor
    %P = diag(sum(S*A.^2,1));
    P = diag(p);
    %B = At*SA;
    %condPB = cond(inv(P)*B)
    [h, flags, relres, iter] = pcg(@(x) (((S*(A*x))')*A)', -g, 1e-6, length(b), P, [], zeros(size(x)));
    iters += iter;

    rr = 1;
    xn = x + h*rr;
    % HACKHACK
    while min(b - A*xn) <= 0
      rr /= 2;
      xn = x + h*rr;
    endwhile
    x = xn;
    ncur = norm(h)/norm(x);
    printf("gnorm=%g ncur=%g\n", gnorm, ncur);
    if ncur < 1e-6 || gnorm < 1e-10
      return
    endif
    k += 1;
  endwhile
endfunction

