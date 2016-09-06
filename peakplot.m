clear all;
figure(1);
exp=load('data');
plot(exp(:,2),exp(:,3),'k');
hold on;
plot(exp(:,2),exp(:,5),'b');
xlim([-5 40]);
hold off;