
figure(1);
data=load('data');
data2=load('TDTRdat.dat');
loglog(data(:,2),data(:,5),'ok');
hold on;
loglog(data2(:,1),data2(:,4),'r','LineWidth',1.5);
hold off;


set(gca, 'FontSize', 14);
set(gca, 'XLim', [100 4000]);
%set(gca, 'YLim', [3 13]);
xlabel('Delay time (ps)');
ylabel('Ratio');
set(gca,'fontsize',14,'linewidth',2);
legend('expt','calc');    