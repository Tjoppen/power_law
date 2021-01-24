% LP solver ála Renegar
if !exist('A','var')
  % load data and reference solution
  printf("loading program.mat\n");
  load('program.mat');
  printf("loading program.m\n");
  program;
  %clear Ai Aj Av;
  program_ref;
  printf("done loading\n");

  A = At';
  clear At;

  % normalize c
  nc = norm(c);
  c = c / nc;

  % augment A with c
  A = [-c'; A];
  b = [0; b];

  % normalize rows
  printf("norming\n");
  d = norm(A, 2, 'rows');
  %d = zeros(size(A,1),1);
  %for k = 1:size(A,1)
  %  d(k) = 1/norm(A(k,:));
  %endfor
  A = diag(d)*A;
  b = diag(d)*b;
  printf("normalized\n");

  % add constraints for xi >= 0
  A = [A; -speye(n)];
  b = [b; zeros(n,1)];

  % compute initial point by going from origin in the direction
  % of c as far as possible. let x0 be the halfway point
  w = w_to_border(A, b, c, zeros(n,1));
  x = w*c*0.5;
  h = d;
  At = A';
  %A2 = A.^2;

  % first get to the initial center
  printf("Computing initial center\n");
  %x = goto_center(A, At, A2, b, x);
  x = goto_center(A, At, b, x);

  nsteps = round(2*129*sqrt(m));
  deltamin = 1/(13*sqrt(m));
  %deltamin = 0.03;
  deltamax = 1/3;
  %delta = 0.03;
  delta = deltamin;
  stats = zeros(nsteps,5);
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

whos
t = time();
tlast = 0;
olast = nc*c'*x/ref;
oblast = olast;

if true
  for step = step0:nsteps
    xbefore = x;
    [x, iters] = goto_center(A, At, b, x);

    %bnext = (1-delta)*b(1) - delta*c'*x;
    %db = bnext - b(1);
    %b(1) = bnext;
    b(1) = (1-delta)*b(1) - delta*c'*x;
    db2 = c'*x + b(1);


    [xnew, iter] = goto_center(A, At, b, x);
    iters += iter;
    h = xnew - x;

    wc = w_to_border(A, b, c, xnew);
    wh = w_to_border(A, b, h, xnew);
    wm = w_to_border(A, b, c+h, xnew);

    a = c'*h/norm(c)/norm(h);
    l = 0.5;

    % solver makes very little progress when cos(theta) approaches 90°
    % put this limit to 0 to make the cutoff at 90°
    amin = 0.05;
    if a > amin
      xh = xnew + 0.5*wh*h;
      xc = xnew + 0.5*wc*c;
      xm = xnew + 0.5*wm*(c+h);
      os = [c'*xh, c'*xc, c'*xm];
      [w, iw] = max(os);
      if iw == 1
        x = xh;
      elseif iw == 2
        x = xc;
      else
        x = xm;
      endif
    else
      x = xnew;
    endif
  
    b(1) = db2 - c'*x;

    tnow = time()-t;
    o = nc*c'*x/ref;
    ob = nc*c'*xbefore/ref;
    %printf("%6i/%6i %.1fsec dt=%.1fsec %g %g %g %f %f %f\n", step, nsteps, tnow, tnow-tlast, nc*c'*x/ref, norm(h), delta, a, wc, wh);
    stats(step,:) = [o, norm(h), delta, iters, tnow-tlast];
    printf("%6i/%6i t=%f dt=%f o=%g (%+g%%) ob=%g (%+g%%) n=%g a=%g i=%i\n", step, nsteps, sum(stats(:,5)), tnow-tlast, o, 100*(o/olast-1), ob, 100*(ob/oblast-1), norm(x-xbefore)/norm(x), a, iters);
    tlast = tnow;
    olast = o;
    oblast = ob;

    if a <= amin
      nsteps = step;
      stats = stats(1:nsteps,:);
      break
    endif
  endfor
else
for step = step0:nsteps
  xbefore = x;
  s = b - A*x;
  if min(s) <= 0
    min(s)
    error("min(s) <= 0");
  endif

  s = 1./s;
  g = At*s;
  S = diag(s)^2;
  SA = S*A;

  if mod(step,100) == 1
    printf("Recalc P\n");
    P = diag(sum(S*A2,1));
    printf("done\n");
  endif

  if k == 1
    x0 = d;
  else
    x0 = zeros(size(x));
  end

  [h, flags, relres, iters] = pcg(@(x) At*(SA*x), -g, 1e-3, 100, P, [], x0);

  x += h;

  dnext = x - xbefore;
  a = d'*dnext/norm(d)/norm(dnext);
  d = dnext;

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
  delta = max(min(delta, deltamax), deltamin);

  tnow = time()-t;
  stats(step,:) = [nc*c'*x/ref, norm(h), delta, iters, tnow-tlast];
  printf("%6i/%6i %.1fsec dt=%.1fsec %g %g %g %i\n", step, nsteps, sum(stats(:,5)), tnow-tlast, nc*c'*x/ref, norm(h), delta, iters);
  tlast = tnow;
  b(1) = (1-delta)*b(1) - delta*c'*x;

  if mod(step,100) == 0
    %semilogy(stats(1:step,3))
    %plot([10*stats(1:step,1),log10(stats(1:step,[2,3]))]);
    %axis([1,step,-5,10*max(stats(1:step,1))])
    %legend('10o','log_{10} |h|','log_{10} \delta');
  endif
endfor
endif
t = time() - t
totiters = sum(stats(:,4))

% m       n       totiters  t       result
perfstats = [
  300,    100,    41045,      25.6, 0;
  1000,   300,    271126,    154.3, 0;
  3000,   1000,   932282,    830.7, 0; % rtol=1e2^k
  3000,   1000,   620797,    538.9, 0; % rtol=1e6
  3000,   1000,   579967,    526.2, 0; % rtol=1e3^k
  3000,   1000,   462974,    436.0, 0; % rtol=1e3
  3000,   1000,   392764,    396.8, 0; % rtol=1e3, x0
  3000,   1000,   392199,    394.8, 0; % norm>30
  3000,   1000,   174539,    209.3, 0.810665; % calm, step mod 20
  3000,   1000,   205079,    240.1, 0.804636; % sqrt(calm)
  3000,   1000,   127988,    176.9, 0.802375; % mod 40
  3000,   1000,   83664,     140.5, 0.823210; % mod 300
  3000,   1000,   81197,     145.4, 0.834342; % mod 1000
  10000,  3000,   77611,     417.3, 0.715042;
  30000,  10000,  279615,   3563.8, 0.628547;
  30000,  10000,  24438,      48.9, 0.977340; % new solver
  100000, 30000,  113705,  19286.4, 0.376237;
  100000, 30000,   47120,   1427.1, 0.974217; % new solver
  300000, 100000, 737021, 117066.2, 0.676991; % 32.5 timmar
  300000, 100000,  20191,   2768.5, 0.902674; % new solver @ step 24 (a=0.256923)
  300000, 100000,  31326,   3623.8, 0.917009; % new solver stop @ step 26 (a=0.0321687)
  1000000,300000,  28425,   6536.8, 0.875182; % herrmann, new solver, stop @ 28 (a=0.0793525), ref guessed
  3000000,1000000, 28891,   2394.6, 0.734454; % laptop, stop @ 30 (a=0.103162), ref guessed
  %7426/141312 1.68623e+13 på 6514.6 sec
  %1891/141312 2629.9sec 1.60163e+13
  % 5148.4sec
];

%plot(stats*diag([1,1e-2,1,1e-3]));

