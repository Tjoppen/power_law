% asdf
d = csvread('all.csv');
d = d(3:8,:);
c = [ones(6,1),log(d(:,1))] \ log(d(:,2));
k = c(1);
m = c(2);
s = diff(log(d(:,1:2)));
mk = max(s(:,2)./s(:,1));

graphics_toolkit('gnuplot')

figure(1)
loglog(d(:,1),d(:,2),'*b')
hold on
grid on
loglog(d(:,1),exp(k+m*log(d(:,1))),'-r')
xlabel('Number of industries')
ylabel('Wall time (s)')
%legend('results',sprintf('m=%f, k=%f', m, k), 'location', 'east')
title(sprintf('Solver wall time measurements and linear fit (slope=%.2f)', m))
print('results-fit.png')

% ./power_law_gen2 200 30 10 50 10 0.5
load('../program.mat')
msz = 4;

figure(2)
spy(A, msz)
title('spy(A)')
print('spy.png')

figure(3)
spy([A;speye(n);c'], msz)
title("spy([A;speye(n);c'])")
print('spy2.png')

figure(4)
spy([A';speye(m);b'], msz)
title("spy([A';speye(m);b'])")
print('spy3.png')
