function [dat] = metac_int_hgf_binary_mab(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')





%% ________________________________________________________________________
%%%% 1 mu_2(0) for all bandits only %%%%


%% fit interoceptive HGF_mab2 modified (Bayes Optimal)

hgf_mab2_config = tapas_hgf_binary_mab_metac_config;
hgf_mab2_config.n_bandits = 2;

for n = 1:size(dat.u_bin,2)

    %take first y_pred as priormu
    hgf_mab2_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1); 
    hgf_mab2_config = tapas_align_priors_fields(hgf_mab2_config)

    dat.hgf_binary_mab.sub(n).bo_est = tapas_fitModel([],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_config,...
        tapas_bayes_optimal_binary_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_plotTraj(dat.hgf_binary_mab.sub(n).bo_est)
    hold on;
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab', ['mab2_bo_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% _______________
% sim interoceptive HGF_mab2 modified (use BO pars)
for n = 1:12
    dat.hgf_binary_mab.sub(n).sim = tapas_simModel([dat.u_bin(:,n) dat.u_mab2],...
        'tapas_hgf_binary_mab_metac',...
        dat.hgf_binary_mab.sub(n).bo_est.p_prc.p,... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
        'tapas_beta_obs',...
        log(128),...
        12345 ...
        );
    tapas_hgf_binary_mab_plotTraj(dat.hgf_binary_mab.sub(n).sim)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab.sub(n).sim.y, '.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab', ['mab2_sim_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% fit interoceptive HGF_mab2 (modified mu_2(0))

hgf_mab2_config = tapas_hgf_binary_mab_metac_config;
hgf_mab2_config.n_bandits = 2;

for n = 1:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_config.priormus = dat.hgf_binary_mab.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_config = tapas_align_priors_fields(hgf_mab2_config)

    dat.hgf_binary_mab.sub(n).est = tapas_fitModel(dat.y_pred(:,n),...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_config,...
        tapas_beta_obs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_metac_plotTraj(dat.hgf_binary_mab.sub(n).est)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab.sub(n).est.y,'.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab', ['mab2_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


end