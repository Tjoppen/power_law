% asdf
function [x, iters] = goto_center (A, b, x0)
% 10k x 3k
%before     totiters =  52795 tottime =  85.736  step 34
%simple     totiters =  61309 tottime =  126.27  step 35
%memory=1   totiters =  51666 tottime =  90.727  step 34

% 30k x 10k
% P         totiters =  101288 tottime =  454.86  step 35
% Pfun      totiters =  67841  tottime =  326.87  step 34
% o=0.95 cutoff:
% Pfun 0.5  totiters =  13628  tottime =  75.181  step 23 o=0.955716 a=0.217663
% Pfun 0.85 totiters =  12119  tottime =  62.749  step 13 o=0.951147 n=0.0510442
% Pfun 0.9  totiters =  11809  tottime =  65.152  step 18 o=0.951329 a=0.238843
% Pfun 0.99 totiters =  19260  tottime =  96.007  step 11 o=0.960487 a=0.228261

  x = x0;
  k = 1;
  iters = 0;
  ss = [];
  ys = [];
  rho = [];

  function y = pr2 (x, P, ss, ys, rho, k)
    if k == 0
      y = P\x;
    else
      z1 = x - ys(:,k)*(rho(k)*(ss(:,k)'*x));
      %z2 = pr2(z1, P, ss, ys, rho, k-1);
      z2 = P\z1;
      z3 = z2 - ss(:,k)*(rho(k)*(ys(:,k)'*z2));
      y = z3 + ss(:,k)*(rho(k)*(ss(:,k)'*x));
    endif
  endfunction

  function y = pr (x, P, ss, ys, rho)
    y = pr2(x, P, ss, ys, rho, size(ss,2));
    %r = x;
    %as = [];
    %for k = 1:size(ss,2)
    %  a = rho(k)*(ss(:,k)'*r);
    %  r -= a*ys(:,k);
    %  as = [as, a];
    %endfor
    %r = P\r;
    %for k = size(ss,2):-1:1
    %  b = rho(k)*(ys(:,k)'*r);
    %  r += (as(k) - b)*ss(:,k);
    %endfor
    %y = r;
  endfunction

  while true
    s = b - A*x;
    if min(s) <= 0
      [w,iw] = min(s)
      k
      error("min(s) <= 0");
    endif
    s = 1./s;
    g = (s'*A)';

    if k > 1
      %ys = [ys, g - glast];
      %rho = [rho, 1/(ys(:,end)'*ss(:,end))];
      ys = g - glast;
      rho = [1/(ys'*ss)];
    endif


    gnorm = norm(g);
    s2 = s.^2;
    S = diag(s2);
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
    elseif flags > 1 % accept it not reaching the specified tolerance
      flags
      error("solve failed");
    endif

    rr = 1;
    xn = x + h*rr;
    while min(b - A*xn) <= 0
      rr /= 2;
      xn = x + h*rr;
    endwhile
    x = xn;
    ncur = norm(h)/norm(x);
    printf("iter=%5i gnorm=%g ncur=%g, cond≃%g\n", iter, gnorm, ncur, max(p)/min(p));
    if ncur < 1e-6 || gnorm < 1e-10
      return
    endif

    %ss = [ss, h];
    ss = h;
    glast = g;

    k += 1;
  endwhile
endfunction

