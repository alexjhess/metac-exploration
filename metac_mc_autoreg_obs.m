function [dat] = metac_mc_autoreg_obs(dat)


%% plot settings
% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% create model space
[mod] = metac_create_model_space(1);


%% fit mc obs models using rawPE as input

dat.F = NaN(size(mod,2),size(dat.u_bin,2));

for m = 1:size(mod,2)

    % init obs param_mat (size n_obs_pars x N)
    mod(m).param_mat = NaN(size(mod(m).obs_idx,2),...
            size(dat.u_pe,2));

    for n = 1:size(dat.u_pe,2)
        
        % fit
        mod(m).sub(n).est = tapas_fitModel(dat.pdat.y_mc(:,n),...
            dat.pdat.u_pe(:,n),...
            mod(m).prc_config,...
            mod(m).obs_config,...
            tapas_quasinewton_optim_config ...
            );

        dat.F(m,n) = mod(m).sub(n).est.optim.LME;
        mod(m).param_mat(:,n) = mod(m).sub(n).est.p_obs.p(mod(m).obs_idx);

        % plot
        figure;
        subplot(2,1,1)
        plot(mod(m).sub(n).est.u.^2)
        ylim([0 1])
        ylabel('raw PE squared')
        subplot(2,1,2)
        if m < 5
            plot(mod(m).sub(n).est.y, '.')
            ylim([0 1])
        else
            plot(log(mod(m).sub(n).est.y ./ (1-mod(m).sub(n).est.y)), '.')
        end
        hold on;
        plot(mod(m).sub(n).est.optim.yhat)
        ylabel('mc response')
        figdir = fullfile('figures', 'mc_autoreg_obs', ['mod' num2str(m) '_est_sub' num2str(n)]);
        print(figdir, '-dpng');
        close;
    end

end


%% model comparison

% collect F matrix
for m = 1:size(mod,2)
    if m < 5
        mod(m) = dat.mod(m)
    end
    for n = 1:size(dat.u_pe,2)
        dat.F(m,n) = mod(m).sub(n).est.optim.LME;
    end
end


%% FFX BMS
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


%% plot residuals
for m = 1:size(mod,2)
    figure;
    for n=1:size(dat.u_pe,2)
        subplot(5,4,n)
        histogram(mod(m).sub(n).est.optim.res)
    end
end


%% extract est params of winning obs model

% npars x N
dat.param_mat = mod(dat.ffx.idx).param_mat;

% plot histogram of est obs param vals
figure
for i = 1:size(dat.param_mat,1)
    subplot(3,2,i)
    histogram(dat.param_mat(i,:))
end


%% save data
close all;

dat.mod = mod;

save(fullfile('data', 'discovery_set_fits_tmp.mat'), 'dat', '-mat');



end




