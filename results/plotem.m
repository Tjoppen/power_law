d=csvread('all.csv');
hold off; loglog(d(:,1),d(:,2),'*'); hold on; grid on; loglog(d(:,1),exp(k+m*log(d(:,1))),'-r')
xlabel('Number of industries'); ylabel('time (s)')
legend('results',sprintf('m=%f, k=%f', m, k), 'location', 'east')
