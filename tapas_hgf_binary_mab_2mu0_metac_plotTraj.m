function tapas_hgf_binary_mab_2mu0_metac_plotTraj(r)
% Plots the estimated or generated trajectories for the binary HGF perceptual model for multi-armed
% bandit situations.
%
% Usage example:  est = tapas_fitModel(responses, inputs); tapas_hgf_binary_plotTraj(est);
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Optional plotting of standard deviations (true or false)
plotsd = true;

% Optional plotting of responses (true or false)
ploty = true;

% Set up display
scrsz = get(0,'screenSize');
outerpos = [0.2*scrsz(3),0.2*scrsz(4),0.8*scrsz(3),0.8*scrsz(4)];
figure(...
    'OuterPosition', outerpos,...
    'Name', 'HGF trajectories');

% Set up colors
colors = [1 0 0; 0.67 0 1; 0 0.67 1; 0.67 1 0];

% Number of bandits
b = r.c_prc.n_bandits;

% Number of trials
n = size(r.u,1);

% Time axis
% if size(r.u,2) > 1
%     t = r.u(:,end)';
% else
    t = ones(1,n);
% end

ts = cumsum(t);
ts = [0, ts];

% Number of levels
try
    l = r.c_prc.n_levels;
catch
    l = (length(r.p_prc.p)+1)/5;
end

% number of plots
if size(r.y,2) <= 1
    npl = l;
else % multimodal obs
    npl = l + size(r.y,2)-1;
end

% Upper levels
for j = 1:l-2

    % Subplots
    subplot(npl,1,j);

    if plotsd == true
        upperprior = r.p_prc.mu_0(l-j+1) +sqrt(r.p_prc.sa_0(l-j+1));
        lowerprior = r.p_prc.mu_0(l-j+1) -sqrt(r.p_prc.sa_0(l-j+1));
        upper = [upperprior; r.traj.mu(:,l-j+1,1)+sqrt(r.traj.sa(:,l-j+1,1))];
        lower = [lowerprior; r.traj.mu(:,l-j+1,1)-sqrt(r.traj.sa(:,l-j+1,1))];
    
        plot(0, upperprior, 'ob', 'LineWidth', 1);
        hold all;
        plot(0, lowerprior, 'ob', 'LineWidth', 1);
        fill([ts, fliplr(ts)], [(upper)', fliplr((lower)')], ...
             'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    end
    plot(ts, [r.p_prc.mu_0(l-j+1); r.traj.mu(:,l-j+1,1)], 'b', 'LineWidth', 2);
    hold all;
    plot(0, r.p_prc.mu_0(l-j+1), 'ob', 'LineWidth', 2); % prior
    xlim([0 ts(end)]);
    title(['Posterior expectation of x_' num2str(l-j+1)], 'FontWeight', 'bold');
    ylabel(['\mu_', num2str(l-j+1)]);
end

% Second level
subplot(npl,1,l-1)


if plotsd == true
    for j=1:b
        if j == 1
            upperprior = r.p_prc.mu_0(2) +sqrt(r.p_prc.sa_0(2));
            lowerprior = r.p_prc.mu_0(2) -sqrt(r.p_prc.sa_0(2));
        elseif j == 2
            upperprior = r.p_prc.mu2_0 +sqrt(r.p_prc.sa_0(2));
            lowerprior = r.p_prc.mu2_0 -sqrt(r.p_prc.sa_0(2));
        end
        upper = [upperprior; r.traj.mu(:,2,j)+sqrt(r.traj.sa(:,2,j))];
        lower = [lowerprior; r.traj.mu(:,2,j)-sqrt(r.traj.sa(:,2,j))];
    
        plot(0, upperprior, 'o', 'Color', colors(j,:), 'LineWidth', 1);
        hold all;
        plot(0, lowerprior, 'o', 'Color', colors(j,:), 'LineWidth', 1);
        fill([ts, fliplr(ts)], [(upper)', fliplr((lower)')], ...
             colors(j,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    end
end
for j=1:b
    plot(ts, [r.p_prc.mu_0(1); r.traj.mu(:,2,j)], 'Color', colors(j,:), 'LineWidth', 2);
    hold all;
    plot(0, r.p_prc.mu_0(1), 'o', 'Color', colors(j,:), 'LineWidth', 2); % prior
end
xlim([0 ts(end)]);
title('Posterior expectations of x_2', 'FontWeight', 'bold');
ylabel('\mu_2');

% Input level
subplot(npl,1,l);

for j=1:b
    if j == 1
        plot(ts, [tapas_sgm(r.p_prc.mu_0(2), 1); tapas_sgm(r.traj.mu(:,2,j), 1)], 'Color', colors(j,:), 'LineWidth', 2);
        hold all;
        plot(0, tapas_sgm(r.p_prc.mu_0(2), 1), 'o', 'Color', colors(j,:), 'LineWidth', 2); % prior
        plot(ts(r.u(:,2)==j)+1, r.u(r.u(:,2)==j,1), '.', 'Color', colors(j,:)); % inputs
    elseif j == 2
        plot(ts, [tapas_sgm(r.p_prc.mu2_0, 1); tapas_sgm(r.traj.mu(:,2,j), 1)], 'Color', colors(j,:), 'LineWidth', 2);
        hold all;
        plot(0, tapas_sgm(r.p_prc.mu2_0, 1), 'o', 'Color', colors(j,:), 'LineWidth', 2); % prior
        plot(ts(r.u(:,2)==j)+1, r.u(r.u(:,2)==j,1), '.', 'Color', colors(j,:)); % inputs
    end
end
% plot(ts(2:end), r.u(:,1), '.', 'Color', [0 0 0]); % inputs
plot(ts(2:end), r.traj.wt(:,1), 'k') % implied learning rate 
if (ploty == true) && ~isempty(find(strcmp(fieldnames(r),'y'))) && ~isempty(r.y)
    y = r.y(:,1);
    if ~isempty(find(strcmp(fieldnames(r),'irr')))
        y(r.irr) = NaN; % weed out irregular responses
        plot(ts(r.irr)+1,  1.08.*ones([1 length(r.irr)]), 'x', 'Color', [1 0.7 0], 'Markersize', 11, 'LineWidth', 2); % irregular responses
    end
    for j=1:b
        % plot(find(y==j), 1.08*ones([1 length(find(y==j))]), '.', 'Color', colors(j,:)); % binary responses
        % plot(ts(r.u(:,2)==j)+1, r.y(r.u(:,2)==j,1), '.', 'Color', colors(j,:)) % responses (color by bandit)
    end
    plot(r.y(:,1), '.', 'Color', [1 0.7 0]) % responses
    title(['Response y, input u (black dots), learning rate (fine black line), and posterior expectation of reward s(\mu_2) ', ...
           '(colour coded) for \rho=', num2str(r.p_prc.rho(2:end)), ', \kappa=', ...
           num2str(r.p_prc.ka(2:end)), ', \omega=', num2str(r.p_prc.om(2:end))], ...
      'FontWeight', 'bold');
    ylabel('y, u, s(\mu_2)');
    axis([0 ts(end) -0.15 1.15]);
else
    title(['Input u (dots), learning rate (fine black line), and posterior expectation of input s(\mu_2) ', ...
           '(red) for \rho=', num2str(r.p_prc.rho(2:end)), ', \kappa=', ...
           num2str(r.p_prc.ka(2:end)), ', \omega=', num2str(r.p_prc.om(2:end))], ...
      'FontWeight', 'bold');
    ylabel('u, s(\mu_2)');
    axis([0 ts(end) -0.1 1.1]);
end
plot(ts(2:end), 0.5, 'k');
xlabel('Trial number');
hold off;

% second response data modality
if npl > 3
    subplot(npl,1,npl);
    plot(r.y(:,2), '.', 'Color', [1 0.7 0]) % mc responses
    hold all;
    if isfield(r, 'c_sim')
        title(['MC response y (orange) for obs param vector ', num2str(r.p_obs.p)], ...
          'FontWeight', 'bold');
    else
        plot(r.optim.yhat(:,2), 'Color', 'b') % mc responses
        
        title(['y_{mc} (orange) and post. expectation of ', ...
               'y_{hat,mc} (blue) for obs transp param vector ', num2str(r.p_obs.ptrans)], ...
          'FontWeight', 'bold');
    end
    xlim([0 n])
    ylabel('$y_{mc}, \hat{y}_{mc}$', 'Interpreter', 'Latex')
    xlabel('Trial number');
    hold off;
end
