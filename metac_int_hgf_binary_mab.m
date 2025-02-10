function [dat] = metac_int_hgf_binary_mab(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% TODO:
%   - SIMULATING FROM MAB models --> adapt OBS FCTs to simulate from
%   different bandits...




%% ________________________________________________________________________
%%%% multiple mu_2(0) (one for each bandit) %%%%

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

windtrials = find(dat.u_mab2==2);

for n = 1:size(dat.u_bin,2)

    %take first y_pred as priormu
    hgf_mab2_2mu0_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1);
    hgf_mab2_2mu0_config.priormus(15) = tapas_logit(dat.y_pred(windtrials(1),n),1);
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config)

    dat.hgf_binary_mab_2mu0.sub(n).bo_est = tapas_fitModel([],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        tapas_bayes_optimal_binary_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0.sub(n).bo_est)
    hold on;
    ts = cumsum(dat.u_mab2');
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0', ['mab2_bo_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% _______________
% sim interoceptive HGF_mab2 with 2 different mu2(0) (use BO pars)

% TODO:
%   - ADAPT OBS FCT TO SIMULATE FROM DIFFERENT BANDITS!

for n = 1%:size(dat.u_bin,2)
    dat.hgf_binary_mab_2mu0.sub(n).sim = tapas_simModel([dat.u_bin(:,n) dat.u_mab2],...
        'tapas_hgf_binary_mab_2mu0_metac',...
        dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.p,... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
        'tapas_beta_obs',...
        log(128),...
        12345 ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0.sub(n).sim)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0.sub(n).sim.y, '.')
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

%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

for n = 1:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_2mu0_config.priormus = dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config)

    dat.hgf_binary_mab_2mu0.sub(n).est = tapas_fitModel(dat.y_pred(:,n),...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        tapas_beta_obs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0.sub(n).est)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0.sub(n).est.y,'.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0', ['mab2_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


%% plot inferred PE (HGF_mab2_2mu0 vs raw PE)

for n = 1:size(dat.u_bin,2)
    figure
    plot(dat.u_pe(:,n))
    hold on;
    plot(dat.hgf_binary_mab_2mu0.sub(n).est.traj.da(:,1))
    % adjust for sgm transform
    adj = dat.hgf_binary_mab_2mu0.sub(n).est.traj.muhat(:,1) .* (1-dat.hgf_binary_mab_2mu0.sub(n).est.traj.muhat(:,1));
    pwPE = dat.hgf_binary_mab_2mu0.sub(n).est.traj.epsi(:,2);
    adj_pwPE = adj .* pwPE;
    % plot(adj_pwPE)
    % color code different conditions
    ax = axis;
    fill([dat.trials, fliplr(dat.trials)],...
        [ax(3)*ones(1,length(dat.trials)), ax(4)*ones(1,length(dat.trials))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    legend('PE_{raw}', 'PE_{HGF mab2 2mu0}')
    figdir = fullfile('figures', ['raw_PE_vs_HGF_mab2_2mu0_PE']);
    print(figdir, '-dpng');
end




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
    tapas_hgf_binary_mab_metac_plotTraj(dat.hgf_binary_mab.sub(n).bo_est)
    hold on;
    ts = cumsum(dat.u_mab2');
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
for n = 1:size(dat.u_bin,2)
    dat.hgf_binary_mab.sub(n).sim = tapas_simModel([dat.u_bin(:,n) dat.u_mab2],...
        'tapas_hgf_binary_mab_metac',...
        dat.hgf_binary_mab.sub(n).bo_est.p_prc.p,... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
        'tapas_beta_obs',...
        log(128),...
        12345 ...
        );
    tapas_hgf_binary_mab_metac_plotTraj(dat.hgf_binary_mab.sub(n).sim)
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

%% fit interoceptive HGF_mab2 (bo_pars including modified mu_2(0))

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