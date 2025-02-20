function [dat] = metac_fit_behaviour(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

%% create model space
[mod, bo] = metac_create_model_space();


%% ________________________________________________________________________
%%%% HGF MAB multiple mu_2(0) (one for each bandit) %%%%
%%%% est Bayes Optimal pars (prc model)

windtrials = find(dat.u_mab2==2);

for n = 1:size(dat.u_bin,2)

    %take first y_pred as priormu
    po.prc_config.priormus(2) = tapas_logit(dat.y_pred(1,n),1);
    bo.prc_config.priormus(15) = tapas_logit(dat.y_pred(windtrials(1),n),1);
    bo.prc_config = tapas_align_priors_fields(bo.prc_config);
    
    % fit
    bo.sub(n).bo_est = tapas_fitModel([],...
        [dat.u_bin(:,n) dat.u_mab2],...
        bo.prc_config,...
        bo.obs_config,...
        tapas_quasinewton_optim_config ...
        );

    % plot
    tapas_hgf_binary_mab_2mu0_metac_plotTraj(bo.sub(n).bo_est)
    figdir = fullfile('figures', 'SBF', ['prc_bo_est_sub' num2str(n)]);
    print(figdir, '-dpng');
    close;

end

%% fit interoceptive HGF_mab2 with different mu2(0) per bandit 
% (bo_pars including modified mu_2(0))
% MC null OBS, MC RES OBS, MC PE OBS, MC RES + PE OBS !

dat.F = NaN(size(mod,2),size(dat.u_bin,2));

for m = 1:size(mod,2)

    % init obs param_mat (size n_obs_pars-1 x N)
    mod(m).param_mat = NaN(size(mod(m).obs_idx(1:end-1),2),...
            size(dat.u_bin,2));

    for n = 1:size(dat.u_bin,2)
    
        %take bo pars as priormus
        mod(m).prc_config.priormus = bo.sub(n).bo_est.p_prc.ptrans; 
        mod(m).prc_config = tapas_align_priors_fields(mod(m).prc_config);
        
        % fit
        mod(m).sub(n).est = tapas_fitModel([dat.y_pred(:,n), dat.y_mc(:,n)],...
            [dat.u_bin(:,n) dat.u_mab2],...
            mod(m).prc_config,...
            mod(m).obs_config,...
            tapas_quasinewton_optim_config ...
            );

        dat.F(m,n) = mod(m).sub(n).est.optim.LME;
        mod(m).param_mat(:,n) = mod(m).sub(n).est.p_obs.p(mod(m).obs_idx(1:end-1));

        % plot
        tapas_hgf_binary_mab_2mu0_metac_plotTraj(mod(m).sub(n).est)
        figdir = fullfile('figures', 'SBF', ['mod' num2str(m) '_est_sub' num2str(n)]);
        print(figdir, '-dpng');
        close;
    end
end


%% model comparison

% FFX BMS
[dat.ffx.sumLME, dat.ffx.pp, dat.ffx.GBF, dat.ffx.ABF] = metac_FFXBMS(dat.F');
[dat.ffx.val, dat.ffx.idx] = max(dat.ffx.sumLME);
disp('FFX BMS results: ')
fprintf('winning model %i \n', dat.ffx.idx)
if ~isempty(dat.ffx.GBF)
    fprintf('GBF %i \n', dat.ffx.GBF)
end
if ~isempty(dat.ffx.ABF)
    fprintf('ABF %i \n', dat.ffx.ABF)
end

%% RFX BMS
[dat.rfx.posterior, dat.rfx.out] = VBA_groupBMC(dat.F);
[dat.rfx.val, dat.rfx.idx] = max(dat.rfx.out.pxp);
disp('RFX BMS results: ')
fprintf('winning model %i \n', dat.rfx.idx)
fprintf('PXP %.2f \n', dat.rfx.out.pxp(dat.rfx.idx))
fprintf('Ef %.2f \n', dat.rfx.out.Ef(dat.rfx.idx))


%% extract est params of winning model

% npars x N
dat.param_mat = mod(dat.ffx.idx).param_mat;


%% save data

dat.mod = mod;

save(fullfile('data', 'discovery_set_fits_tmp.mat'), 'dat', '-mat');



end