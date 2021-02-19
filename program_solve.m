% LP solver experiment
if !exist('A','var')
  [A, b, c, nc, m, n, ref] = program_load('program.mat');
  x = program_init(A, b, c);

  nsteps = round(2*129*sqrt(m));
  delta = 1/(13*sqrt(m));
  stats = zeros(nsteps,5);
endif

if exist('step','var')
  printf("Continuing @ step %i\n", step);
  step0 = step;
else
  printf("Performing %i steps with delta=%g\n", nsteps, delta);
  step0 = 1;
endif

olast = nc*c'*x/ref;
oblast = olast;

for step = step0:nsteps
  t = time();
  xbefore = x;
  [x, iters] = goto_center(A, b, x);

  b(1) = (1-delta)*b(1) - delta*c'*x;
  db2 = c'*x + b(1);


  [xnew, iter] = goto_center(A, b, x);
  iters += iter;
  h = xnew - x;
  a = c'*h/norm(c)/norm(h);

  % solver makes very little progress when cos(theta) approaches 90°
  % put this limit to 0 to make the cutoff at 90°
  amin = 0.00;
  if a > amin
    bestl = 1;
    bests = -Inf;
    bestw = 0;
    printf("Trying different directions ");
    nh = h/norm(h);
    Ax = A*xnew;
    for l = exp(-0.005:0.00025:0.005)
      printf(".");
      fflush(stdout);
      hc = l*nh + (1-l)*c;
      w = w_to_border(A, b, hc, Ax);
      s = w*(c'*hc);
      if s > bests
        bests = s;
        bestl = l;
        bestw = w;
        besthc = hc;
      endif
    endfor
    printf(" l=%f was best\n", bestl);
    % go forward quite aggressively
    x = xnew + 0.85*bestw*besthc;
    b(1) = db2 - c'*x;
  else
    x = xnew;
    b(1) = (1-delta)*b(1) - delta*c'*x;
  endif

  dt = time() - t;
  o = nc*c'*x/ref;
  ob = nc*c'*xbefore/ref;
  stats(step,:) = [o, norm(h), delta, iters, dt];

  printf("%6i/%6i t=%f dt=%f o=%g (%+g%%) ob=%g (%+g%%) n=%g a=%g i=%i\n", step, nsteps, sum(stats(:,5)), dt, o, 100*(o/olast-1), ob, 100*(ob/oblast-1), norm(x-xbefore)/norm(x), a, iters);

  olast = o;
  oblast = ob;
  [res, xx, oo] = check_solution(A(2:end,:), b(2:end), c, x);

  if res % || a < min
    printf("Found optimal solution!\n");
    x = xx;
    stats(step,:) = [nc*oo/ref, norm(h), delta, iters, dt];
    nsteps = step;
    stats = stats(1:nsteps,:);
    break
  endif
endfor

totiters = sum(stats(:,4))
tottime = sum(stats(:,5))

% m       n       totiters     t    o/ref
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
  300000, 100000, 737021, 117066.2, 0.676991; % 32.5 hours
  300000, 100000,  20191,   2768.5, 0.902674; % new solver @ step 24 (a=0.256923)
  300000, 100000,  31326,   3623.8, 0.917009; % new solver @ step 26 (a=0.0321687)
  1000000,300000,  28425,   6536.8, 0.875182; % herrmann, new solver, stop @ 28 (a=0.0793525), ref guessed
  3000000,1000000, 28891,   2394.6, 0.734454; % laptop, stop @ 30 (a=0.103162), ref guessed
  10000000,3000000,14109,  38838.4, 0;
];

