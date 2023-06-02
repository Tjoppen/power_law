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
function w = w_to_border_dual (A, c, dl, ll, b, od)
  epsilon = 1e-6;
  num = c - (ll'*A)';
  den = -(dl'*A)';
  %li = find(abs(den) > epsilon);
  num = [num; od - b'*ll];
  den = [den; b'*dl];
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

  l = -ll ./ dl;
  l = find(l > epsilon);
  w = min([w; l]);
endfunction
