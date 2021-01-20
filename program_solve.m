% load data and reference solution
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
d = [];
for k = 1:size(A,1)
  d = [d; 1/norm(A(k,:))];
  %A(k,:) /= nn;
  %b(k) /= nn;
  %[k,m]
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
%h = x*0.1;
w
norm(x)

function x = conjgrad(A, S, b, x)
    r = b - A' * (S * (A * x));
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A' * (S * (A * p));
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
    i
end

d = zeros(size(x));
h = d;
At = A';
kmax = round((10*sqrt(m)));

for k = 1:kmax
  xx = x;
  while true
    s = b - A*x;
    if min(s) <= 0
      min(s)
      error("min(s) <= 0");
    endif
    g = A'*(1./s);
    S = inv(diag(s))^2;
    SA = S*A;
    if true %k <= 3
      B = At*SA;
      h = B\-g;
    else
      %h = conjgrad(A, S, -g, d); %zeros(size(h)));
      %h = bicgstab(@(x) At*(SA*x), -g, 1e-3, 10*length(g), [], [], h);
      h = pcg(@(x) At*(SA*x), -g, 1e-6, 10*length(g), [], [], d*0.5);
      %h = cgs(@(x) At*(SA*x), -g, 1e-6, 10*length(g), [], [], d);
    endif

    %[min(s), max(s), norm(s), norm(x), norm(g), norm(h), norm(h)/norm(x)]
    rr = 1;
    xn = x + h*rr;
    % HACKHACK
    while min(b - A*xn) <= 0
      rr /= 2;
      xn = x + h*rr;
    endwhile
    x = xn;
    norm(h*rr)/norm(x)
    if norm(h*rr)/norm(x) < 1e-6
      break
    endif
  endwhile
  d = x - xx;
  [100*k/kmax, nc*c'*x/ref, -b(1), nc*c'*x]
  b(1) -= (c'*x + b(1))*0.33;
  %b(1) -= (c'*x + b(1))/10;
endfor
[nc*c'*x/ref, -b(1)]
