function [dat] = metac_int_hgf_binary_mc_autoreg_obs_mab(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% TODO:
%   - SIMULATING FROM MAB models --> adapt OBS FCTs to simulate from
%   different bandits...




%% ________________________________________________________________________
%%%% HGF MAB multiple mu_2(0) (one for each bandit) %%%%
%% est Bayes Optimal pars (prc model)

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
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_bo_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% _______________
% sim interoceptive HGF_mab2 with 2 different mu2(0) (use BO pars)
% MC null OBS model

% TODO:
%   - ADAPT OBS FCT TO SIMULATE FROM DIFFERENT BANDITS!

for n = 1%:size(dat.u_bin,2)
    dat.hgf_binary_mab_2mu0_mc_null.sub(n).sim = tapas_simModel([dat.u_bin(:,n) dat.u_mab2],...
        'tapas_hgf_binary_mab_2mu0_metac',...
        dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.p,... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
        'mc_null_inobs',...
        [0.05 0.5 0.005],...
        12345 ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_null.sub(n).sim)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0_mc_null.sub(n).sim.y, '.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_null_sim_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC null OBS !

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

for n = 1:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_2mu0_config.priormus = dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config);

    dat.hgf_binary_mab_2mu0_mc_null.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_null_inobs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_null.sub(n).est)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0_mc_null.sub(n).est.y,'.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_null_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC RES OBS !

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

for n = 1:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_2mu0_config.priormus = dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config);

    dat.hgf_binary_mab_2mu0_mc_res.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_autoreg_res_inobs_config,... mc_autoreg_pe_inobs_config
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_res.sub(n).est)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0_mc_res.sub(n).est.y,'.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_res_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC PE OBS !

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

for n = 1%:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_2mu0_config.priormus = dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config);

    dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_autoreg_pe_inobs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est)
    hold on;
    ts = cumsum(dat.u_mab2');
    plot(ts,dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est.y,'.')
    % color code conditions
    ax=axis;
    fill([ts, fliplr(ts)],...
        [ax(3)*ones(1,length(ts)), ax(4)*ones(1,length(ts))],...
        [dat.u_mab4', fliplr(dat.u_mab4')], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    colormap(flipud(colormap("autumn")))
    hold off;
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_pe_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end









%% ________________________________________________________________________
% build autoregr mc obs model (giuliara). HGF MAB ???
n=1;
dat.ehgf_binary_mc_pe.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_autoreg_pe_inobs',... 'mc_null_inobs',...
    [0.05 1 -2 -1.5 -0.04 0.005],... [0.05 -2 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_pe.sub(n).sim)
hold on
plot(dat.ehgf_binary_mc_pe.sub(n).sim.y(:,2), '.')

%% sim null model. HGF MAB ???
dat.ehgf_binary_mc_null.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_null_inobs',...
    [0.05 0.5 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_null.sub(n).sim)
hold on
plot(dat.ehgf_binary_mc_null.sub(n).sim.y(:,2), '.')


%% fit metacognitive HGF (giuliara). HGF MAB ???

for n = 1:size(dat.u_bin,2)
    dat.ehgf_binary_mc_null.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        dat.u_bin(:,n),...
        tapas_ehgf_binary_config,...
        mc_null_inobs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_null.sub(n).est)
    hold on
    plot(dat.ehgf_binary_mc_null.sub(n).est.y(:,2), '.')
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ____________
% build autoregr mc obs model (giuliara). NOT HGF MAB 
n=1;
dat.ehgf_binary_mc_pe.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_autoreg_pe_inobs',... 'mc_null_inobs',...
    [0.05 1 -2 -1.5 -0.04 0.005],... [0.05 -2 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_pe.sub(n).sim)
hold on
plot(dat.ehgf_binary_mc_pe.sub(n).sim.y(:,2), '.')

%% sim null model. NOT HGF MAB 
dat.ehgf_binary_mc_null.sub(n).sim = tapas_simModel(dat.u_bin(:,n),...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_null_inobs',...
    [0.05 0.5 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_null.sub(n).sim)
hold on
plot(dat.ehgf_binary_mc_null.sub(n).sim.y(:,2), '.')


%% fit metacognitive HGF (giuliara). NOT HGF MAB 

for n = 1%:size(dat.u_bin,2)
    dat.ehgf_binary_mc_null.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        dat.u_bin(:,n),...
        tapas_ehgf_binary_config,...
        mc_null_inobs_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_null.sub(n).est)
    hold on
    plot(dat.ehgf_binary_mc_null.sub(n).est.y(:,2), '.')
end

%% bo pars

for n = 1%:size(dat.u_bin,2)
    dat.ehgf_binary_mc_null.sub(n).est = tapas_fitModel([],...
        dat.u_bin(:,n),...
        tapas_ehgf_binary_config,...
        tapas_bayes_optimal_binary_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_ehgf_binary_plotTraj(dat.ehgf_binary_mc_null.sub(n).est)
end

