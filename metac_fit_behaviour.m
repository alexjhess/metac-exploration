function [dat] = metac_fit_behaviour(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% ________________________________________________________________________
%%%% HGF MAB multiple mu_2(0) (one for each bandit) %%%%
%%%% est Bayes Optimal pars (prc model)

% 2 different mu2(0)
hgf_mab2_2mu0_config = tapas_hgf_binary_mab_2mu0_metac_config;
hgf_mab2_2mu0_config.n_bandits = 2;

windtrials = find(dat.u_mab2==2);

for n = 1:size(dat.u_bin,2)

    %take first y_pred as priormu
    hgf_mab2_2mu0_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1);
    hgf_mab2_2mu0_config.priormus(15) = tapas_logit(dat.y_pred(windtrials(1),n),1);
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config);

    dat.hgf_binary_mab_2mu0.sub(n).bo_est = tapas_fitModel([],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        tapas_bayes_optimal_binary_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0.sub(n).bo_est)
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_bo_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC null OBS !

for n = 1:size(dat.u_bin,2)

    %take bo pars as priormus
    hgf_mab2_2mu0_config.priormus = dat.hgf_binary_mab_2mu0.sub(n).bo_est.p_prc.ptrans; 
    hgf_mab2_2mu0_config = tapas_align_priors_fields(hgf_mab2_2mu0_config);

    dat.hgf_binary_mab_2mu0_mc_null.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_null_inobs_mab_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_null.sub(n).est)
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_null_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end

%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC RES OBS !

for n = 1:size(dat.u_bin,2)
    dat.hgf_binary_mab_2mu0_mc_res.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_autoreg_res_inobs_config,... mc_autoreg_pe_inobs_config
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_res.sub(n).est)
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_res_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end



%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC PE OBS !

for n = 1:size(dat.u_bin,2)
    dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_autoreg_pe_inobs_mab_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est)
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_pe_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


%% fit interoceptive HGF_mab2 with different mu2(0) per bandit (bo_pars including modified mu_2(0))
% MC PE+RES OBS !

for n = 1:size(dat.u_bin,2)
    dat.hgf_binary_mab_2mu0_mc_pe_res.sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
        [dat.u_bin(:,n) dat.u_mab2],...
        hgf_mab2_2mu0_config,...
        mc_autoreg_pe_res_inobs_mab_config,...
        tapas_quasinewton_optim_config ...
        );
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(dat.hgf_binary_mab_2mu0_mc_pe_res.sub(n).est)
    figdir = fullfile('figures', 'int_hgf_binary_mab_2mu0_mc_autoreg', ['mab2_mc_pe_res_est_sub' num2str(n)]);
    print(figdir, '-dpng');
end


%% model comparison

dat.F = NaN(4,size(dat.u_bin,2));

for n = 1:size(dat.u_bin,2)
    dat.F(1,n) = dat.hgf_binary_mab_2mu0_mc_null.sub(n).est.optim.LME;
    dat.F(2,n) = dat.hgf_binary_mab_2mu0_mc_pe.sub(n).est.optim.LME;
    dat.F(3,n) = dat.hgf_binary_mab_2mu0_mc_res.sub(n).est.optim.LME;
    dat.F(4,n) = dat.hgf_binary_mab_2mu0_mc_pe_res.sub(n).est.optim.LME;
end

% FFX BMS
[dat.ffx.sumLME, dat.ffx.pp, dat.ffx.GBF] = metac_FFXBMS(dat.F');
[dat.ffx.val, dat.ffx.idx] = max(dat.ffx.sumLME);
disp('FFX BMS results: ')
fprintf('winning model %i \n', dat.ffx.idx)
if ~isempty(dat.ffx.GBF)
    fprintf('GBF %i \n', dat.ffx.GBF)
end

% RFX BMS
[dat.rfx.posterior, dat.rfx.out] = VBA_groupBMC(dat.F);
[dat.rfx.val, dat.rfx.idx] = max(dat.rfx.out.pxp);
disp('RFX BMS results: ')
disp(sprintf('winning model %i', dat.rfx.idx))
fprintf('PXP %.2f \n', dat.rfx.out.pxp(dat.rfx.idx))
fprintf('Ef %.2f \n', dat.rfx.out.Ef(dat.rfx.idx))


%% extract est params of winning model

% TODO: insert winning model !!!!!!

% npars x N
dat.param_mat = NaN(size(dat.hgf_binary_mab_2mu0_mc_pe_res.sub(n).est.p_obs.p,2),...
    size(dat.u_bin,2));

for n = 1:size(dat.u_bin,2)
    dat.param_mat(:,n) = dat.hgf_binary_mab_2mu0_mc_pe_res.sub(n).est.p_obs.p;
end


%% centralize values.. ?

param_mat_std = NaN(size(dat.param_mat));


%% save data
save(fullfile('data', 'discovery_set_fits_tmp.mat'), 'dat', '-mat');



end