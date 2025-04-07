function [dat] = metac_init_priors_grid_search(dat, pdat, mod, n_pe)

%% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

% seed for rng
rng(123, 'twister')
options.rng.settings = rng;
options.rng.idx = 1; % Set counter for random number states

% create model space
% [mod, bo] = metac_create_model_space(1); % 1=logit space

% n_pe = 1; % PE traj from pilot sub n_pe
%%
for m=1:size(mod,2)
    %%
    input.prc.transInp = mod(m).prc_config.priormus;
    c.c_prc = mod(m).prc_config;
    input.prc.nativeInp = mod(m).prc_config.transp_prc_fun(c, input.prc.transInp);
    c.c_obs = mod(m).obs_config;
    obs_vect_nat = c.c_obs.transp_obs_fun(c, c.c_obs.priormus);
    obs_idx = mod(m).obs_idx;
    
    n_sim = 100;
    free_param = NaN(n_sim,size(mod(m).obs_config.priormus,2));
    
    if m >= 4
        % zemc
        free_param(:,1) = 0:4/99:4;
        % gamma
        free_param(:,2) = -1:2/99:1;
        % Shift
        free_param(:,3) = -3:10/99:7;
        % Scale
        free_param(:,4) = -5:10/99:5;
        % y0
        free_param(:,5) = -1:2/99:1;
        % w
        free_param(:,6) = -2:4/99:2;
    elseif m == 1
        % zemc
        free_param(:,1) = 0:4/99:4;
        % Shift
        free_param(:,2) = -5:10/99:5;
    else
        % zemc
        free_param(:,1) = 0:4/99:4;
        % gamma
        free_param(:,2) = -1:2/99:1;
        % Shift
        free_param(:,3) = -3:10/99:7;
        % Scale
        free_param(:,4) = -5:10/99:5;
        % y0
        free_param(:,5) = -1:2/99:1;
    end
    
    
    
    disp(obs_vect_nat)
    priormu_sim = tapas_simModel([dat.u_bin(:,n_pe) pdat.u_pe(:,n_pe)],...
        mod(m).prc,...
        input.prc.nativeInp,...
        mod(m).obs,...
        obs_vect_nat,...
        options.rng.settings.State(1));
    
    %%
    for k = 1:size(free_param,2)
    %%
        for n=1:n_sim
        
            % get priormus
            input.obs.transInp = mod(m).obs_config.priormus;
            
            % change value of param k
            input.obs.nativeInp = mod(m).obs_config.transp_obs_fun(c, input.obs.transInp);
            input.obs.nativeInp(k) = free_param(n,k); %input.obs.nativeInp(mod(m).obs_idx(k)) = free_param(n,k);
         
            % sim
            sub(n,m,k).sim = tapas_simModel([dat.u_bin(:,n_pe) pdat.u_pe(:,n_pe)],...
                mod(m).prc,...
                input.prc.nativeInp,...
                mod(m).obs,...
                input.obs.nativeInp,...
                options.rng.settings.State(options.rng.idx, 1));
            
            % save sim input
            sim.sub(n,m,k).input = input;
            
            % Update the rng state idx
            options.rng.idx = options.rng.idx+1;
            if options.rng.idx == (length(options.rng.settings.State)+1)
                options.rng.idx = 1;
            end
        
        end

        % reset rng idx
        options.rng.idx = 1;
        
        
        %% plot simulated responses
        figure
        hold on;
        for n = 1:n_sim
            if k == 1
                col = [0 (1-free_param(n,k)/max(free_param(:,k))) 0.5 ...
                    (1-free_param(n,k)/max(free_param(:,k)))];
            elseif k == 2
                col = [0 0.5 0.5 ...
                    0.5];
            else
                col = [0 0.5*(1-free_param(n,k)/max(free_param(:,k))) 0.5 ...
                    0.5*(1-free_param(n,k)/max(free_param(:,k)))];
            end
            plot(sub(n,m,k).sim.y, 'Color', col)
        end
        
        %% u_pe_pos = 0.5*(dat.u_pe + 1.0);
        u_pe_sq = 0.95.*(pdat.u_pe_sq-0.5)+0.5;
        logit_pe_sq = log(u_pe_sq ./ (1-u_pe_sq));
        if m <= 4
            plot(pdat.u_pe_sq(:,n_pe), 'r', 'LineWidth',2)
            % plot(logit_pe_sq(:,n_pe), 'y', 'LineWidth',2)
        elseif m>4
            plot(pdat.u_pe(:,n_pe), 'r', 'LineWidth',2)
        end
    
        pdat.logit_y_mc = log(pdat.y_mc ./ (1-pdat.y_mc));
        % plot(mean(pdat.logit_y_mc,2), 'k', 'LineWidth',2)
        plot(pdat.y_mc(:,n_pe), 'k', 'LineWidth',2)
        plot(priormu_sim.y, 'k', 'LineWidth', 2, 'LineStyle','--')
        hold off;

        figdir = fullfile('figures', 'logit_mc_autoreg_obs', 'gridsearch_priors',...
            ['linesearch_mod' num2str(m) '_nPE', num2str(n_pe), '_param' num2str(k)]);
        print(figdir, '-dpng');
        close;
    
    end

end


end