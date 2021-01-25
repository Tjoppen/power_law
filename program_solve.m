% LP solver ála Renegar
if !exist('A','var')
  % load data and reference solution
  printf("loading program.mat\n");
  load('program.mat');
  program_ref;
  printf("done loading\n");

  A = A';
  %clear At;

  % normalize c
  nc = norm(c);
  c = c / nc;

  % augment A with c
  A = [-c'; A];
  b = [0; b];

  % normalize rows
  printf("norming\n");
  d = 1./norm(A, 2, 'rows');
  D = diag(d);
  A = D*A;
  b = D*b;
  clear d D;
  printf("normalized\n");

  % add constraints for xi >= 0
  A = [A; -speye(n)];
  b = [b; zeros(n,1)];

  % compute initial point by going from origin in the direction
  % of c as far as possible. let x0 be the halfway point
  w = w_to_border(A, b, c, zeros(size(A,1),1));
  x = w*c*0.5;
  %h = d;
  %At = A';
  %A2 = A.^2;

  % first get to the initial center
  printf("Computing initial center\n");
  %x = goto_center(A, At, A2, b, x);
  x = goto_center(A, b, x);

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
  whos
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
olast = nc*c'*x/ref;
oblast = olast;

if true
  for step = step0:nsteps
    xbefore = x;
    [x, iters] = goto_center(A, b, x);

    %bnext = (1-delta)*b(1) - delta*c'*x;
    %db = bnext - b(1);
    %b(1) = bnext;
    b(1) = (1-delta)*b(1) - delta*c'*x;
    db2 = c'*x + b(1);


    [xnew, iter] = goto_center(A, b, x);
    iters += iter;
    h = xnew - x;
    a = c'*h/norm(c)/norm(h);

    % solver makes very little progress when cos(theta) approaches 90°
    % put this limit to 0 to make the cutoff at 90°
    amin = 0.05;
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
      x = xnew + 0.5*bestw*besthc;
    else
      x = xnew;
    endif
  
    b(1) = db2 - c'*x;

    tnow = time()-t;
    o = nc*c'*x/ref;
    ob = nc*c'*xbefore/ref;
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
%15/815868 t=20885.029821 dt=1712.381353 o=0.381574 (+8.1391%) ob=0.352855 (+6.42727%) n=0.150753 a=0.278966 i=392
%16/815868 t=22865.568467 dt=1980.538646 o=0.414083 (+8.51981%) ob=0.381574 (+8.1391%) n=0.141528 a=0.293602 i=457
%17/815868 t=25189.984377 dt=2324.415910 o=0.451483 (+9.03181%) ob=0.414083 (+8.51981%) n=0.135194 a=0.310758 i=542
%18/815868 t=27899.674790 dt=2709.690413 o=0.493395 (+9.28336%) ob=0.451483 (+0%) n=0.127614 a=0.328527 i=634
%19/815868 t=31756.931224 dt=3857.256434 o=0.536083 (+8.65178%) ob=0.493395 (+0%) n=0.112574 a=0.344733 i=760
%25/815868 t=66737.991659 dt=7367.324870 o=0.730429 (-3.74017e-06%) ob=0.730429 (+4.73781%) n=0.00564382 a=0.0425988 i=1520
%t =  38838.40404
%totiters =  14109

];

%plot(stats*diag([1,1e-2,1,1e-3]));

