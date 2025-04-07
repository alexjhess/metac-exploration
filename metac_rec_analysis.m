function [rec] = metac_rec_analysis(vs, rec, n_sim)

%% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

%% model identifiability (LME Winner classification)

% pre-allocate
class.LMEwinner = NaN(size(vs.mod, 2), size(vs.mod, 2));
class.percLMEwinner = NaN(size(class.LMEwinner));

% calc winner freq for each data generating model
for m = 1:size(vs.mod, 2)
    [class.max(m).val, class.max(m).idx] = max(rec.model(m).LME, [], 2);
    for i = 1:size(vs.mod, 2)
        class.LMEwinner(m,i) = sum(class.max(m).idx==i);
    end
    class.percLMEwinner(m,:) = class.LMEwinner(m,:)./n_sim;
    % accuracy
    class.acc(m) = class.percLMEwinner(m,m);
end

% balanced accuraccy
class.balacc = mean(class.acc);
% chance threshold (inv binomial distr)
class.chancethr = binoinv(0.9, n_sim, 1/size(vs.mod, 2)) / n_sim;
% save to struct
rec.class = class;

%% parameter recovery (Pearson's correlation coefficient)
for m = 1:size(vs.mod, 2)
    % % prc model
    % [prc_coef, prc_p] = corr(rec.param.prc(m).sim, rec.param.prc(m).est);
    % rec.param.prc(m).pcc = diag(prc_coef);
    % rec.param.prc(m).pval = diag(prc_p);
    % obs model
    [obs_coef, obs_p] = corr(rec.param.obs(m).sim, rec.param.obs(m).est);
    rec.param.obs(m).pcc = diag(obs_coef);
    rec.param.obs(m).pval = diag(obs_p);
end

%% model identifiability (RFX BMS)

% pre-allocate
bmc.rfx.Ef = NaN(size(vs.mod, 2), size(vs.mod, 2));
bmc.rfx.ep = NaN(size(bmc.rfx.Ef));
bmc.rfx.pxp = NaN(size(bmc.rfx.Ef));

% toolbox settings for BMS results display
bmc.opt.verbose = false;
bmc.opt.DisplayWin = false;

% run BMS for each data generating model
for m = 1:size(vs.mod, 2)
    L = rec.model(m).LME';
    [bmc.post, bmc.out] = VBA_groupBMC(L, bmc.opt); %evtl add options...
    bmc.rfx.Ef(m,:) = bmc.out.Ef';
    bmc.rfx.ep(m,:) = bmc.out.ep;
    bmc.rfx.pxp(m,:) = bmc.out.pxp;
end

% save to struct
rec.bmc = bmc;

%% FIG: plot param rec
for m = 1:size(vs.mod,2)
    npars = length(vs.mod(m).prc_idx)+length(vs.mod(m).obs_idx);
    nx = ceil(sqrt(npars));
    ny = round(sqrt(npars));
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    for j = 1:npars
        subplot(nx, ny, j)
        if j > length(vs.mod(m).prc_idx) %% all obs param
            k = j-length(vs.mod(m).prc_idx);
            scatter(rec.param.obs(m).sim(:,k), rec.param.obs(m).est(:,k), 15, 'k', 'filled');
            pcc = rec.param.obs(m).pcc(k);
            hline = refline(1,0);
            hline.Color = 'k';
            % if k == 1
            %     title('log(\zeta)')
            % elseif k == size(vs.mod(m).obs_idx, 2)
            %     title('log(\Sigma)')
            % elseif k == 2
            %     title('\beta_0')
            % elseif k == 3
            %     title('\beta_1')
            % elseif k == 4
            %     title('\beta_2')
            % elseif k == 5
            %     title('\beta_3')
            % elseif k == 6
            %     title('\beta_4')
            % end
        else
            % scatter(rec.param.prc(m).sim(:,j), rec.param.prc(m).est(:,j), 15, 'k', 'filled');
            % pcc = rec.param.prc(m).pcc(j);
            % hline = refline(1,0);
            % hline.Color = 'k';
            % if j == 1
            %     title('\omega_2')
            % elseif j == 2
            %     title('\omega_3')
            % end
        end        
        xlabel('Simulated values');
        ylabel('Recovered values');
        str = sprintf('r = %1.2f', pcc);
        textXpos = min(get(gca, 'xlim')) + (max(get(gca, 'xlim')) - min(get(gca, 'xlim')))*0.05;
        textYpos = max(get(gca, 'ylim')) - (max(get(gca, 'ylim')) - min(get(gca, 'ylim')))*0.05;
        T = text(textXpos, textYpos, str);
        set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    end
    figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
        ['param_rec_m', num2str(m)]);
        print(figdir, '-dpng');
        close;
end


%% FIG: model ident (LME winner classification)

% axis labels
for m = 1:size(vs.mod,2)
    mod_nr{m} = ['M', num2str(m)];
end

figure
bwr = @(n)interp1([1 2], [[1 1 1]; [0 0 0]], linspace(1, 2, n), 'linear');
imagesc(rec.class.percLMEwinner)
colormap(bwr(64));
colorbar;
set(gca, 'clim', [0 1])
ax = gca; ax.FontSize = 20;
ax.XTick = [1:size(vs.mod,2)];
ax.XTickLabel = mod_nr;
ax.YTick = [1:size(vs.mod,2)];
ax.YTickLabel = mod_nr;

figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
    ['model_rec_LME']);
print(figdir, '-dpng');
close;

%% FIG: model ident (rfx bms, PXP)

% axis labels
for m = 1:size(vs.mod,2)
    mod_nr{m} = ['M', num2str(m)];
end

figure
bwr = @(n)interp1([1 2], [[1 1 1]; [0 0 0]], linspace(1, 2, n), 'linear');
imagesc(rec.bmc.rfx.pxp)
colormap(bwr(64));
colorbar;
set(gca, 'clim', [0 1])
ax = gca; ax.FontSize = 20;
ax.XTick = [1:size(vs.mod,2)];
ax.XTickLabel = mod_nr;
ax.YTick = [1:size(vs.mod,2)];
ax.YTickLabel = mod_nr;

figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'ep',...
    ['model_rec_PXP']);
print(figdir, '-dpng');
close;



end