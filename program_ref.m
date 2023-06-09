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
% computed by lp_solve
%ref = 56061.19888294;

%ref = 1434158.69711017; % m=1000, n=300, real    0m14,342s
% 100.00000         0.95007       452.98517   1362554.31329
% 22.940 seconds
% 100.00000         0.95010       452.99660   1362588.56194
% 198.06 s med tol=1e-10
% 100.00000         0.94995       452.92499   1362373.87638
% 36.588 s med tol=1e-3

%ref = 1898028.53073738; % m=3000, n=1000, real    30m8,327s

%ref = 2.18685925; % ./a.out 10 3 10
%ref = 56061.19888294;  % ./a.out 300 100 10
%ref = 1434158.69711017; % ./a.out 1000 300 10
%ref = 56412635.63686652; % ./a.out 3000 1000 10
ref = 1398990377.15945506; % ./a.out 10000 3000 10
%ref = 53102724588.38802338; % ./a.out 30000 10000 10
%ref = 1415019154359.47387695; % time ./a.out 100000 30000 10 = 1656.78user
%ref = 53648789919600.83593750; % time ./a.out 300000 100000 10 = 44930.27user
%ref = 1;
%ref = 1415019154359.47387695*1010.9; % 1M 300k guess
%ref = 53648789919600.83593750*1010.9; % 3M 1M guess
%ref = 4.60921583204459e+16; % 3M 1M result (not lp_solve)
%ref = 1.05622586096742e+18; % 10M 30M result (not lp_solve)

