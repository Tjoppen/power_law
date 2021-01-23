% LP solver ála Renegar
if !exist('A','var')
  % load data and reference solution
  printf("loading program\n");
  program;
  program_ref;
  printf("loaded\n");

  % normalize c
  nc = norm(c);
  c = c / nc;

  % augment A with c
  A = [-c'; A];
  b = [0; b];

  % normalize rows
  printf("norming\n");
  d = zeros(size(A,1),1);
  for k = 1:size(A,1)
    d(k) = 1/norm(A(k,:));
  endfor
  A = diag(d)*A;
  b = diag(d)*b;
  printf("normalized\n");

  % add constraints for xi >= 0
  A = [A; -speye(n)];
  b = [b; zeros(n,1)];

  epsilon = 1e-6;

  % compute initial point by going from origin in the direction
  % of c as far as possible. let x0 be the halfway point
  num = b; % - A*x;
  den = A*c;
  %li = find(abs(den) > epsilon);
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

  x = w*c*0.5;
  h = d;
  At = A';
  A2 = A.^2;

  % first get to the initial center
  printf("Computing initial center\n");
  nprev = 0;
  rprev = 1;
  while true
    s = b - A*x;
    if min(s) <= 0
      min(s)
      error("min(s) <= 0");
    endif
    g = At*(1./s);
    S = inv(diag(s))^2;
    SA = S*A;
    %B = At*SA;
    %h = B\-g;
    P = diag(sum(S*A2,1));
    rtol = 1e6;
    [h, flags, relres, iter] = pcg(@(x) At*(SA*x), -g, norm(d)/rtol, 10*length(g), P, [], zeros(size(x)));

    rr = 1;
    xn = x + h*rr;
    % HACKHACK
    while min(b - A*xn) <= 0
      rr /= 2;
      xn = x + h*rr;
    endwhile
    x = xn;
    ncur = norm(h*rr)/norm(x);
    ncur
    rcur = nprev / ncur;
    if ncur < 1e-6 %&& rcur > 30*rprev
      break
    endif
    nprev = ncur;
    rprev = rcur;
  endwhile

  nsteps = round(2*129*sqrt(m));
  %deltamin = 1/(13*sqrt(m));
  deltamin = 0.03;
  deltamax = 1/3;
  delta = 0.03;
  stats = zeros(nsteps,4);
  calm = 0;
  vroom = 0;
  d = x*0.5;
endif

if exist('step','var')
  printf("Continuing @ step %i\n", step);
  step0 = step;
else
  printf("Performing %i steps with delta=%g\n", nsteps, delta);
  step0 = 1;
endif

t = time();
tlast = 0;
for step = step0:nsteps
  hs = [];
  nk = 2;

  if calm >= 3
    if vroom
      printf("vroom!\n");
      vroom = 0;
    endif
    nk = 1;
    if mod(step,50) == 1
      nk = 2;
    endif
  end

  xbefore = x;
  iters = 0;
  for k = 1:nk
    s = b - A*x;
    if min(s) <= 0
      min(s)
      error("min(s) <= 0");
    endif

    s = 1./s;
    g = At*s;
    S = diag(s)^2;
    SA = S*A;

    if mod(step,100) == 1 && k == 1
      printf("Recalc P\n");
      P = diag(sum(S*A2,1));
      printf("done\n");
    endif

    if k == 1
      x0 = d;
    else
      x0 = zeros(size(x));
    end

    [h, flags, relres, iter] = pcg(@(x) At*(SA*x), -g, 1e-3, 100, P, [], x0);

    iters += iter;
    x += h;
    hs = [hs; norm(h)];
  endfor
  dnext = x - xbefore;
  a = d'*dnext/norm(d)/norm(dnext);
  d = dnext;

  % TODO: justera delta beroende på antal iterationer i pcg()
  if step > 10 && nk > 1
    if hs(1)/norm(h) < 10
      %delta *= 0.9;
      %calm = 0;
      %printf("back off!\n");
    else
      if calm == 0
        printf("calm..\n");
        vroom = 1;
      endif
      calm += 1;
      if hs(1)/norm(h) > 30
        delta *= 1.01;
      endif
    endif
  elseif nk == 1
    if iters < 2
      delta *= 1.02;
    elseif iters < 5
      delta *= 1.01;
    elseif iters > 30
      delta *= 0.8;
    elseif iters > 20
      delta *= 0.9;
    elseif iters > 15
      delta *= 0.99;
    %elseif iters > 8
    %  delta *= 0.99;
    end
  endif
  delta = max(min(delta, deltamax), deltamin);

  tnow = time()-t;
  printf("%6i/%6i %.1fsec dt=%.1fsec %g %g %g %i\n", step, nsteps, tnow, tnow-tlast, nc*c'*x/ref, hs(1), delta, iters);
  tlast = tnow;
  stats(step,:) = [nc*c'*x/ref, hs(1), delta, iters];
  b(1) = (1-delta)*b(1) - delta*c'*x;

  if mod(step,100) == 0
    %semilogy(stats(1:step,3))
    %plot([10*stats(1:step,1),log10(stats(1:step,[2,3]))]);
    %axis([1,step,-5,10*max(stats(1:step,1))])
    %legend('10o','log_{10} |h|','log_{10} \delta');
  endif
endfor
t = time() - t
totiters = sum(stats(:,4))

% m      n     totiters  t       result
perfstats = [
  300,   100,  41045,    25.578, 0;
  1000,  300,  271126,   154.36, 0;
  3000,  1000, 932282,   830.71, 0; % rtol=1e2^k
  3000,  1000, 620797,   538.90, 0; % rtol=1e6
  3000,  1000, 579967,   526.29, 0; % rtol=1e3^k
  3000,  1000, 462974,   436.06, 0; % rtol=1e3
  3000,  1000, 392764,   396.86, 0; % rtol=1e3, x0
  3000,  1000, 392199,   394.80, 0; % norm>30
  3000,  1000, 174539,   209.30, 0.810665; % calm, step mod 20
  3000,  1000, 205079,   240.10, 0.804636; % sqrt(calm)
  3000,  1000, 127988,   176.93, 0.802375; % mod 40
  3000,  1000, 83664,    140.51, 0.823210; % mod 300
  3000,  1000, 81197,    145.48, 0.834342; % mod 1000
  10000, 3000, 77611,    417.35, 0.715042;
  30000, 10000,279615,   3563.8, 0.628547;
  100000,30000,113705,  19286.4, 0.376237;
  300000,100000,737021,117066.2, 0.676991; % 32.5 timmar
  %7426/141312 1.68623e+13 på 6514.6 sec
  %1891/141312 2629.9sec 1.60163e+13
  % 5148.4sec
];

plot(stats*diag([1,1e-2,1,1e-3]));

