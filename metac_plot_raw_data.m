function [] = metac_plot_raw_data(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

%% plot task variables (sub 1)
n=1;
figure
ax1 = plot(dat.u_bin(:,n), '.');
hold on;
plot(dat.y_pred(:,n))
plot(dat.y_mc(:,n))
ylim([-0.1 1.1])
% color code different phases
ax = axis;
fill([dat.trials, fliplr(dat.trials)],...
    [ax(3)*ones(1,length(dat.trials)), ax(4)*ones(1,length(dat.trials))],...
    [dat.phase, fliplr(dat.phase)], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
colormap(flipud(colormap("autumn")))
legend({'u', '$y_{int}$', '$y_{mc}$', 'phase'},...
    'Interpreter','latex', 'Position',[0.75 0.25 0.15 0.1])
hold off;
figdir = fullfile('figures', ['task_vars_phase_sub', num2str(n)]);
print(figdir, '-dpng');

%% plot task variables (sub 1)
n=1;
figure
ax1 = plot(dat.u_bin(:,n), '.');
hold on;
plot(dat.y_pred(:,n))
plot(dat.y_mc(:,n))
ylim([-0.1 1.1])
% color code different conditions
ax = axis;
fill([dat.trials, fliplr(dat.trials)],...
    [ax(3)*ones(1,length(dat.trials)), ax(4)*ones(1,length(dat.trials))],...
    [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
colormap(flipud(colormap("autumn")))
legend({'u', '$y_{int}$', '$y_{mc}$', 'condition'},...
    'Interpreter','latex', 'Position',[0.75 0.25 0.15 0.1])
hold off;
figdir = fullfile('figures', ['task_vars_cond_sub', num2str(n)]);
print(figdir, '-dpng');

%% plot conditions
figure
ax = plot(dat.u_mab2, '.');
hold on;
plot(dat.u_mab4, 'o');
ylim([0.9 4.1])
legend({'no wind / wind', 'n+s / n+l / w+s / w+l'}, 'Location', 'West')
figdir = fullfile('figures', ['task_cond']);
print(figdir, '-dpng');

%% plot y_mc & avg pos + neg PEs
figure;
plot(mean(dat.y_mc,2), 'LineWidth',2)
hold on;
plot(1-mean(dat.u_pe_pos,2,'omitnan'))
plot(mean(dat.u_pe_neg,2,'omitnan'))
% color code different conditions
ax = axis;
fill([dat.trials, fliplr(dat.trials)],...
    [ax(3)*ones(1,length(dat.trials)), ax(4)*ones(1,length(dat.trials))],...
    [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
colormap(flipud(colormap("autumn")))
hold off;
legend('avg u_{bin}', '1-avg u_{PE+}', 'avg u_{PE-}', 'Location','northeast')
figdir = fullfile('figures', ['avg_control_vs_raw_PE_pos_neg']);
print(figdir, '-dpng');

%% plot avg responses
figure;
plot(mean(dat.y_mc,2))
hold on;
plot(mean(dat.u_pe,2))
% color code different conditions
ax = axis;
fill([dat.trials, fliplr(dat.trials)],...
    [ax(3)*ones(1,length(dat.trials)), ax(4)*ones(1,length(dat.trials))],...
    [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
colormap(flipud(colormap("autumn")))
hold off;
legend('avg y_{mc}', 'avg u_{PE}')
figdir = fullfile('figures', ['avg_control_vs_raw_PE']);
print(figdir, '-dpng');


% %% plot raw data
% figure;
% plot(dat.y_mc)
% hold on;
% plot(1-abs(dat.u_pe))
% plot(dat.y_pred)
% plot(dat.u_bin, '.')
% hold off;
% legend('y_{mc}', '1-|u_{PE}| (raw)', 'y_{int}', 'u_{bin}')
% 
% disp('correlation between y_mc and 1-abs(u_pe):')
% disp(corr(dat.y_mc,1-abs(dat.u_pe)))

% %% plot recoded data
% u_bin_recoded = (dat.u_bin.*2)-1;
% y_pred_recoded = (dat.y_pred*2)-1;
% 
% figure;
% plot(u_bin_recoded, '.')
% hold on;
% plot(y_pred_recoded)
% plot(u_bin_recoded-y_pred_recoded)
% hold off;
% legend('u_{bin}', 'y_{int}', 'PE')


end