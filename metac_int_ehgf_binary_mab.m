function [dat] = metac_int_ehgf_binary_mab(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% fit interoceptive HGF_mab modified (Bayes Optimal)

ehgf_mab_config = tapas_ehgf_ar1_binary_mab_config_metac;
ehgf_mab_config.n_bandits = 4;

for n = 1%:size(dat.u_bin,2)

    %take first y_pred as priormu
    ehgf_mab_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1); 
    ehgf_mab_config = tapas_align_priors_fields(ehgf_mab_config)

    dat.hgf_binary_mab.sub(n).bo_est = tapas_fitModel([],...
        [dat.u_bin(:,n) dat.u_mab4],...
        ehgf_mab_config,...
        tapas_bayes_optimal_binary_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_plotTraj(dat.hgf_binary_mab.sub(n).bo_est)
end

%% _______________
% sim interoceptive HGF_mab modified

for n = 1%:12
    dat.hgf_binary_mab.sub(n).bo_est.p_prc.p
    dat.hgf_binary_mab.sub(n).sim = tapas_simModel([dat.u_bin(:,n) dat.u_mab2],...
        'tapas_hgf_binary_mab_metac',...
        [NaN 2 1 NaN 0.1000 1 NaN 0 0 1 1 NaN -8.1733 -6.0410],... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
        'tapas_beta_obs',...
        log(128),...
        12345 ...
        );
    tapas_hgf_binary_mab_plotTraj(dat.hgf_binary_mab.sub(n).sim)
    hold on;
    plot(dat.hgf_binary_mab.sub(n).sim.y)
    hold off;
end

%% fit interoceptive HGF_mab (modified mu_2(0))

ehgf_mab_config = tapas_ehgf_ar1_binary_mab_config_metac;
ehgf_mab_config.n_bandits = 2;

hgf_mab_config = tapas_hgf_binary_mab_metac_config;
hgf_mab_config.n_bandits = 2;
hgf_mab_config.priormus(13) = -1;
hgf_mab_config.priormus(14) = -2;
hgf_mab_config = tapas_align_priors_fields(hgf_mab_config);

for n = 1%:size(dat.u_bin,2)

    hgf_mab_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1); %take first y_pred as priormu
    hgf_mab_config = tapas_align_priors_fields(hgf_mab_config)

    dat.hgf_binary_mab.sub(n).est = tapas_fitModel(dat.y_pred(:,n),...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab_config,...
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
    figdir = fullfile('figures', 'int_ehgf_binary_mab', ['est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


end