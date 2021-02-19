% Loads the program, normalizes and sets up A
function [A, b, c, nc, m, n, ref] = program_load(name)
  % load data and reference solution
  printf("loading\n");
  load(name);
  program_ref;
  printf("done loading\n");

  % normalize c
  nc = norm(c);
  c = c / nc;

  % un-transpose, augment A with c, add constraints for xi >= 0
  A = [-c'; A'; -speye(n)];
  b = [0; b; zeros(n,1)];

  % normalize rows
  printf("norming\n");
  d = 1./norm(A, 2, 'rows');
  D = diag(d);
  A = D*A;
  b = D*b;
  clear d D;
  printf("normalized\n");
endfunction
