% asdf
d1 = csvread('herrman-results/all.csv');
d2 = csvread('herrman-results-price/all.csv');
n1 = csvread('nnz.csv');
n2 = csvread('price-nnz.csv');
r1 = 2:8;
r2 = [2,3,4,5,6,8];
v1 = d1(r1,1);
v2 = d2(r2,1);
d1(:,1) = n1(:,2);
d2(:,1) = n2(:,2);
d1 = d1(r1,:);
d2 = d2(r2,:);
c1 = [ones(size(d1,1),1),log(d1(:,1))] \ log(d1(:,2));
c2 = [ones(size(d2,1),1),log(d2(:,1))] \ log(d2(:,2));
%k = c(1);
%m = c(2);
%s = diff(log(d(:,1:2)));
%mk = max(s(:,2)./s(:,1));

%graphics_toolkit('gnuplot')

close all
figure(1)
hold on
grid on
loglog(d1(:,1),d1(:,2),'+b')
loglog(d1(:,1),exp(c1(1)+c1(2)*log(d1(:,1))),'-b')
loglog(d2(:,1),d2(:,2),'xr')
loglog(d2(:,1),exp(c2(1)+c2(2)*log(d2(:,1))),'--r')
legend('interconnected',sprintf('%.0f \\mu{}s * nnz(S)^{%.2f}', exp(c1(1))*1e6, c1(2)),'Price',sprintf('%.f \\mu{}s * nnz(S)^{%.2f}', exp(c2(1))*1e6, c2(2)),'location','southeast')
xlabel('nnz(S)')
ylabel('Wall time (s)')
title('Solver wall time measurements and power law fits')
print('results-fit.png')

figure(2)
semilogx(v1, d1(:,3), '-+b', v2, d2(:,3), '--xr')
grid on
axis([1e2,1e6,0,7000]);
legend('interconnected', 'Price','location','southeast')
xlabel('v')
ylabel('Number of CG iterations')
title('Number of CG steps taken vs number of industries')
print('results-steps.png')

figure(3)
semilogx(v1, 100*d1(:,5)./d1(:,4) - 100, '-+b', v2, 100*d2(:,5)./d2(:,4) - 100, '--xr')
grid on
axis([1e2,1e6,0,7]);
legend('interconnected', 'Price','location','southeast')
xlabel('v')
ylabel('Final duality gap (%)')
title('Final duality gap vs number of industries')
print('results-gap.png')

