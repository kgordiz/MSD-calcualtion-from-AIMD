% Linear fit to MSD vs. Deltat data obtained from AIMD simulations

clear

msd_data = importdata('MSD_vs_dt_40_1100C.txt');
x = msd_data(:,1);
y = mean(msd_data(:,2:4),2);

begin_data = 1;
end_data = 2e5; 

figure('color','w')
plot(x/1000,y,'b','linewidth',1.5)

y_for_fit = y(begin_data:end_data);
x_for_fit = x(begin_data:end_data);
p = polyfit(x_for_fit,y_for_fit,1);
p(1) % 10 is to convert from A^2/fs to cm^2/s
hold on
plot([x_for_fit(1),x_for_fit(end)]/1000,[p(1)*x_for_fit(1)+p(2),p(1)*x_for_fit(end)+p(2)],'r--','linewidth',1.5)
set(gca,'fontsize',15)
axis square
xlabel('\Deltat(ps)')
ylabel('MSD (A^2)')
%ylabel('MSD ($\AA^{2}$)','Interpreter','latex','Fontname','helvetica')
