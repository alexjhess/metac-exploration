function [dat] = metac_fit_behaviour2(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

%% create model space
[mod, bo] = metac_create_model_space(1); % 1=logit space


%% fit mc obs models using rawPE as input

dat.main.F = NaN(size(mod,2),size(dat.u_bin,2));
dat.main.Ll = NaN(size(mod,2),size(dat.u_bin,2));
dat.main.comp = NaN(size(mod,2),size(dat.u_bin,2));

for m = 1:size(mod,2)

    % overwrite with pilot priors
    mod(m).obs_config = dat.main.ModSpace(m).obs_config;

    % init obs param_mat (size n_obs_pars x N)
    mod(m).param_mat = NaN(size(mod(m).obs_idx,2),...
            size(dat.u_pe,2));

    for n = 1:size(dat.u_pe,2)
        
        % fit
        mod(m).sub(n).est = tapas_fitModel(dat.pdat.y_mc(:,n),...
            [dat.u_bin(:,n) dat.pdat.u_pe(:,n)],...
            mod(m).prc_config,...
            mod(m).obs_config,...
            tapas_quasinewton_optim_config ...
            );

        dat.main.F(m,n) = mod(m).sub(n).est.optim.LME;
        dat.main.Ll(m,n) = mod(m).sub(n).est.optim.accu;
        dat.main.comp(m,n) = mod(m).sub(n).est.optim.comp;
        mod(m).param_mat(:,n) = mod(m).sub(n).est.p_obs.p(mod(m).obs_idx);

        % plot
        figure;
        subplot(2,1,1)
        if m == 5
            u_mc = mod(m).sub(n).est.u(:,1) + mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            hold on;
            PEsq = 0.95.*(mod(m).sub(n).est.u(:,2).^2-0.5)+0.5; % Shrink 1/2 by a factor of 0.95 (s.t. in *open* unit interval)
            logit_PEsq = log(PEsq ./ (1-PEsq)); % logit transform raw PE squared4
            plot(logit_PEsq)
            legend('R', 'logit PEsq')
            ylabel('logit raw PE squared + Res')
        elseif m == 2 % Res
            u_mc = mod(m).sub(n).est.u(:,1) + mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            ylabel('Res')
            ylim([-2 2])
        elseif m == 3 % PE
            PE = 1+mod(m).sub(n).est.u(:,2); % make sure PE is >0
            logit_PE = log(PE ./ (2-PE)); % logit transform raw PE squared
            plot(logit_PE)
            ylabel('logit raw PE')
        elseif m == 4 % PE + Res
            u_mc = mod(m).sub(n).est.u(:,1) + mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            hold on;
            PE = 1+mod(m).sub(n).est.u(:,2); % make sure PE is >0
            logit_PE = log(PE ./ (2-PE)); % logit transform raw PE squared
            plot(logit_PE)
            legend('R', 'logit PE')
            ylabel('logit raw PE + Res')
        end
        % ylim([0 1])
        subplot(2,1,2)
        plot(log(mod(m).sub(n).est.y ./ (1-mod(m).sub(n).est.y)), '.')
        hold on;
        plot(mod(m).sub(n).est.optim.yhat)
        ylabel('mc response')
        figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
            ['mod' num2str(m) '_est_sub' num2str(n)]);
        print(figdir, '-dpng');
        close;
    end

end


%% model comparison

%% plot acc + comp term of LME
figure
subplot(2,1,1)
bar(dat.Ll')
ylim([-200,0])
title('accuracy')
subplot(2,1,2)
bar(dat.comp')
xlabel('sub')
ylim([0,200])
title('complexity')
legend('null', 'res', 'pe', 'full')
figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep', ['LME_decomp_acc_comp']);
print(figdir, '-dpng');
close;


%% FFX BMS
[dat.main.ffx.sumLME, dat.main.ffx.pp, dat.main.ffx.GBF, dat.main.ffx.ABF] = metac_FFXBMS(dat.main.F');
[dat.main.ffx.val, dat.main.ffx.idx] = max(dat.main.ffx.sumLME);
disp('FFX BMS results: ')
fprintf('winning model %i \n', dat.main.ffx.idx)
if ~isempty(dat.main.ffx.GBF)
    fprintf('GBF %i \n', dat.main.ffx.GBF)
end
if ~isempty(dat.main.ffx.ABF)
    fprintf('ABF %i \n', dat.main.ffx.ABF)
end

%% RFX BMS
[dat.main.rfx.posterior, dat.main.rfx.out] = VBA_groupBMC(dat.main.F);
[dat.main.rfx.val, dat.main.rfx.idx] = max(dat.main.rfx.out.pxp);
disp('RFX BMS results: ')
fprintf('winning model %i \n', dat.main.rfx.idx)
fprintf('PXP %.2f \n', dat.main.rfx.out.pxp(dat.main.rfx.idx))
fprintf('Ef %.2f \n', dat.main.rfx.out.Ef(dat.main.rfx.idx))


%% extract est params of winning model

% npars x N
dat.main.param_mat = mod(dat.main.ffx.idx).param_mat;


%% save data

dat.main.mod = mod;

save(fullfile('data', 'discovery_set_fits_ep_tmp.mat'), 'dat', '-mat');



end