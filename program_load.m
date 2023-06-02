%    power_law
%    Copyright (C) 2023  Tomas Härdin
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as
%    published by the Free Software Foundation, either version 3 of the
%    License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.
%
%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
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
