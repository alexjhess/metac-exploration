function [gest, mod] = metac_fit_behaviour(dat, pdat, mod, priors)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

% create model space
% [mod, bo] = metac_create_model_space(1); % 1=logit spacemod = dat.mod;


%% fit mc obs models using rawPE as input

gest.F = NaN(size(mod,2),size(dat.u_bin,2));
gest.Ll = NaN(size(mod,2),size(dat.u_bin,2));
gest.comp = NaN(size(mod,2),size(dat.u_bin,2));

for m = 1:size(mod,2)

    % init obs param_mat (size n_obs_pars x N)
    mod(m).param_mat = NaN(size(mod(m).obs_idx,2),...
            size(dat.u_pe,2));

    for n = 1:size(dat.u_pe,2)
        
        % fit
        mod(m).sub(n).est = tapas_fitModel(pdat.y_mc(:,n),...
            [dat.u_bin(:,n) pdat.u_pe(:,n)],...
            mod(m).prc_config,...
            mod(m).obs_config,...
            tapas_quasinewton_optim_config ...
            );

        gest.F(m,n) = mod(m).sub(n).est.optim.LME;
        gest.Ll(m,n) = mod(m).sub(n).est.optim.accu;
        gest.comp(m,n) = mod(m).sub(n).est.optim.comp;
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
        if priors == 1
            figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ip',...
                ['mod' num2str(m) '_est_sub' num2str(n)]);
        elseif priors == 2
            figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
                ['mod' num2str(m) '_est_sub' num2str(n)]);
        end
        print(figdir, '-dpng');
        close;
    end

end


%% model comparison

%% plot acc + comp term of LME
figure
subplot(2,1,1)
bar(gest.Ll')
ylim([-200,0])
title('accuracy')
subplot(2,1,2)
bar(gest.comp')
xlabel('sub')
ylim([0,200])
title('complexity')
legend('null', 'res', 'pe', 'full')
if priors == 1
    figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ip',...
        ['LME_decomp_acc_comp']);
elseif priors == 2
    figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
        ['LME_decomp_acc_comp']);
end
print(figdir, '-dpng');
close;


%% FFX BMS
[gest.ffx.sumLME, gest.ffx.pp, gest.ffx.GBF, gest.ffx.ABF] = metac_FFXBMS(gest.F');
[gest.ffx.val, gest.ffx.idx] = max(gest.ffx.sumLME);
disp('FFX BMS results: ')
fprintf('winning model %i \n', gest.ffx.idx)
if ~isempty(gest.ffx.GBF)
    fprintf('GBF %i \n', gest.ffx.GBF)
end
if ~isempty(gest.ffx.ABF)
    fprintf('ABF %i \n', gest.ffx.ABF)
end

%% RFX BMS
[gest.rfx.posterior, gest.rfx.out] = VBA_groupBMC(gest.F);
[gest.rfx.val, gest.rfx.idx] = max(gest.rfx.out.pxp);
disp('RFX BMS results: ')
fprintf('winning model %i \n', gest.rfx.idx)
fprintf('PXP %.2f \n', gest.rfx.out.pxp(gest.rfx.idx))
fprintf('Ef %.2f \n', gest.rfx.out.Ef(gest.rfx.idx))
if priors == 1
    figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ip',...
        ['rfx_bms']);
elseif priors == 2
    figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
        ['rfx_bms']);
end
print(figdir, '-dpng');
close;


%% extract est params of winning model

% npars x N
gest.param_mat = mod(gest.ffx.idx).param_mat;



end