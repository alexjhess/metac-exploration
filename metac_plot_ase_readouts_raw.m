function [] = metac_plot_ase_readouts_raw(dat)

%% avg value per sub
figure
plot(dat.task.avg_c)%, '.')
hold on;
plot(dat.task.avg_tol)%, '.')
plot(dat.task.avg_av)%, '.')
hold off;
legend('control', 'tol', 'av')
ylim([0 1])
xlabel('sub')

%% indiv sub
figure
plot(dat.y_c(10:10:80,:))
hold on;
plot(dat.y_tol(10:10:80,:))
plot(dat.y_av(10:10:80,:))
hold off;
ylim([0 1])

%% mean over sub
figure
avg_c = mean(dat.y_c, 2, 'omitnan');
avg_tol = mean(dat.y_tol, 2, 'omitnan');
avg_av = mean(dat.y_av, 2, 'omitnan');
plot(avg_c(10:10:80))
hold on;
plot(avg_tol(10:10:80))
plot(avg_av(10:10:80))
hold off;
ylim([0 1])
xlabel('trial')
legend('control', 'tol', 'av')

end