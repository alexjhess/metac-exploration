function [] = metac_job_plot(EULER)
% [] = metac_job_plot(EULER)
%
% Creates figures for all steps of the analysis pipeline.
%
% INPUT
%   EULER        binary           Binary indicator variable
%
%   OPTIONAL:
%
% OUTPUT    
%   argout       type
%
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2025 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________

%% setup path
sdir = metac_paths(EULER);

%% load analysis specifications
opts = load(fullfile(sdir, 'results', ['options']));

%% load final results
% res = load(fullfile(sdir, 'results', 'main', ['main_results']));
ds = load(fullfile(sdir, 'results', 'ds', ['discovery_set_fits_tmp.mat']));

%% ________________________________________________________________________
% start producing figures...

%% plot ASE-inspired readouts
metac_plot_ase_readouts_raw(ds.dat);

%% DS: fitted trajs

for m = 1:size(ds.ip.mod,2)
    for n = 1:opts.ds.nS
        
        % plot
        figure;
        subplot(2,1,1)
        if m == 5
            u_mc = ds.mod(m).sub(n).est.u(:,1) + ds.mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            hold on;
            PEsq = 0.95.*(ds.mod(m).sub(n).est.u(:,2).^2-0.5)+0.5; % Shrink 1/2 by a factor of 0.95 (s.t. in *open* unit interval)
            logit_PEsq = log(PEsq ./ (1-PEsq)); % logit transform raw PE squared4
            plot(logit_PEsq)
            legend('R', 'logit PEsq')
            ylabel('logit raw PE squared + Res')
        elseif m == 2 % Res
            u_mc = ds.mod(m).sub(n).est.u(:,1) + ds.mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            ylabel('Res')
            ylim([-2 2])
        elseif m == 3 % PE
            PE = 1+ds.mod(m).sub(n).est.u(:,2); % make sure PE is >0
            logit_PE = log(PE ./ (2-PE)); % logit transform raw PE squared
            plot(logit_PE)
            ylabel('logit raw PE')
        elseif m == 4 % PE + Res
            u_mc = ds.mod(m).sub(n).est.u(:,1) + ds.mod(m).sub(n).est.u(:,1) - 1; % recode to {-1,1}
            plot(u_mc)
            hold on;
            PE = 1+ds.mod(m).sub(n).est.u(:,2); % make sure PE is >0
            logit_PE = log(PE ./ (2-PE)); % logit transform raw PE squared
            plot(logit_PE)
            legend('R', 'logit PE')
            ylabel('logit raw PE + Res')
        end
        % ylim([0 1])
        subplot(2,1,2)
        plot(log(ds.mod(m).sub(n).est.y ./ (1-ds.mod(m).sub(n).est.y)), '.')
        hold on;
        plot(ds.mod(m).sub(n).est.optim.yhat)
        ylabel('mc response')
        % save fig
        figdir = fullfile('figures', 'ds', 'ip',...
            ['mod' num2str(m) '_est_sub' num2str(n)]);
        print(figdir, '-dpng');
        close;
    end

end

%% DS: plot acc + comp term of LME
figure
subplot(2,1,1)
bar(ds.gest.Ll')
ylim([-200,0])
title('accuracy')
subplot(2,1,2)
bar(ds.gest.comp')
xlabel('sub')
ylim([0,200])
title('complexity')
legend('null', 'res', 'pe', 'full')
figdir = fullfile('figures', 'ds', 'ip',...
    ['LME_decomp_acc_comp']);
print(figdir, '-dpng');
close;

%% DS: plot correlations
figure
gplotmatrix(table2array(ds.tab),[],[],[],[],[],[],[],ds.tab.Properties.VariableNames)
figdir = fullfile('figures', 'ds', ['scatterplotmat_ds']);
print(figdir, '-dpng');
close;



%% PILOT: plot estimated priors
% for m = 1:size(res.main.ModSpace,2)
%     npars = length(res.main.ModSpace(m).prc_idx)+length(res.main.ModSpace(m).obs_idx);
%     nx = ceil(sqrt(npars));
%     ny = round(sqrt(npars));
%     figure
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     for j = 1:npars
%         subplot(nx, ny, j)
%         if j > length(res.main.ModSpace(m).prc_idx) % obs
%             k = j - length(res.main.ModSpace(m).prc_idx);
%             idx = res.main.ModSpace(m).obs_idx(k);
%             if k == 1 %log(zeta)
%                 x_min = -3;
%                 x_max = 8;
%             elseif k == size(res.pilot.priors.mod(m).obs_est, 2)
%                 x_min = res.pilot.ModSpace(m).obs_config.priormus(idx)...
%                     -3*res.pilot.ModSpace(m).obs_config.priorsas(idx);
%                 x_max = res.pilot.ModSpace(m).obs_config.priormus(idx)...
%                     +3*res.pilot.ModSpace(m).obs_config.priorsas(idx);
%             else %betas
%                 x_min = -30;
%                 x_max = 30;
%             end
%             x = x_min:0.1:x_max;
% 
%             y = normpdf(x, res.main.ModSpace(m).obs_config.priormus(idx), sqrt(res.main.ModSpace(m).obs_config.priorsas(idx)));
%             y_prior = normpdf(x, res.pilot.ModSpace(m).obs_config.priormus(idx), sqrt(res.pilot.ModSpace(m).obs_config.priorsas(idx)));
% 
%             plot(x, y, 'k')
%             hold on
%             plot(x, y_prior, 'k--')
%             plot(res.pilot.priors.mod(m).obs_est(:,k), -0.05, 'ko')
%             ylim([-0.1 1])
%             str = sprintf('mu = %1.2f, Sa = %1.2f', res.main.ModSpace(m).obs_config.priormus(idx), res.main.ModSpace(m).obs_config.priorsas(idx));
%             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%             set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%             if k == 1
%                 title('log(\zeta)')
%             elseif k == size(res.pilot.priors.mod(m).obs_est, 2)
%                 title('log(\Sigma)')
%             elseif k == 2
%                 title('\beta_0')
%             elseif k == 3
%                 title('\beta_1')
%             elseif k == 4
%                 title('\beta_2')
%             elseif k == 5
%                 title('\beta_3')
%             elseif k == 6
%                 title('\beta_4')
%             end
%         else %prc
%             idx = res.main.ModSpace(m).prc_idx(j);
%             x_min = res.pilot.ModSpace(m).prc_config.priormus(idx)...
%                 -3*res.pilot.ModSpace(m).prc_config.priorsas(idx);
%             x_max = res.pilot.ModSpace(m).prc_config.priormus(idx)...
%                 +3*res.pilot.ModSpace(m).prc_config.priorsas(idx);
%             x = x_min:0.1:x_max;
% 
%             y = normpdf(x, res.main.ModSpace(m).prc_config.priormus(idx), sqrt(res.main.ModSpace(m).prc_config.priorsas(idx)));
%             y_prior = normpdf(x, res.pilot.ModSpace(m).prc_config.priormus(idx), sqrt(res.pilot.ModSpace(m).prc_config.priorsas(idx)));
% 
%             plot(x, y, 'k')
%             hold on
%             plot(x, y_prior, 'k--')
%             plot(res.pilot.priors.mod(m).prc_est(:,j), -0.05, 'ko')
%             ylim([-0.1 1])
%             str = sprintf('mu = %1.2f, Sa = %1.2f', res.main.ModSpace(m).prc_config.priormus(idx), res.main.ModSpace(m).prc_config.priorsas(idx));
%             T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%             set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%             if j == 1
%                 title('\omega_2')
%             elseif j == 2
%                 title('\omega_3')
%             end
%         end
% 
%     end
%     legend('pilot prior', 'initial prior', 'MAP estimates', 'Position', [0.94 0.48 0.03 0.07])
%     figdir = fullfile(sdir, 'figures', 'pilots', 'priors',...
%         ['Priors_model', num2str(m)]);
%     print(figdir, '-dtiff');
%     close;
% end


%% print message
disp('All figures created and saved.')

end