function [sim] = metac_sim_ep(dat, pdat, mod, n_pe, n_sim)

%% dock all figures
set(0,'DefaultFigureWindowStyle','docked')

% seed for rng
rng(123, 'twister')
options.rng.settings = rng;
options.rng.idx = 1; % Set counter for random number states

% pre-allocate variables
sim = struct();
input = struct();

% n_sim = 100;

% loop over synthetic subjects & model space
for n = 1:n_sim
   
    for m = 1:size(mod,2)
        
        % sample free parameter values
        input.prc.transInp = mod(m).prc_config.priormus;
        input.obs.transInp = mod(m).obs_config.priormus;
        for j = 1:size(mod(m).prc_idx,2)
            input.prc.transInp(mod(m).prc_idx(j)) = ...
                normrnd(mod(m).prc_config.priormus(mod(m).prc_idx(j)),...
                abs(sqrt(mod(m).prc_config.priorsas(mod(m).prc_idx(j)))));
        end
        for k = 1:size(mod(m).obs_idx,2)
            input.obs.transInp(mod(m).obs_idx(k)) = ...
                normrnd(mod(m).obs_config.priormus(mod(m).obs_idx(k)),...
                abs(sqrt(mod(m).obs_config.priorsas(mod(m).obs_idx(k)))));
        end
        
        % create simulation input vectors (native space)
        c.c_prc = mod(m).prc_config;
        input.prc.nativeInp = mod(m).prc_config.transp_prc_fun(c, input.prc.transInp);
        c.c_obs = mod(m).obs_config;
        input.obs.nativeInp = mod(m).obs_config.transp_obs_fun(c, input.obs.transInp);
        
        % simulate responses
        % stable = 0;
        % while stable == 0
            % me = [];
            % try
                sim.sub(n,m).data = tapas_simModel([dat.u_bin(:,n_pe) pdat.u_pe(:,n_pe)],...
                    mod(m).prc,...
                    input.prc.nativeInp,...
                    mod(m).obs,...
                    input.obs.nativeInp,...
                    options.rng.settings.State(options.rng.idx, 1));
                stable = 1;
            % catch me
            %     fprintf('simulation failed for Model %1.0f, synth. Sub %1.0f \n', [m, n]);
            %     fprintf('Prc Param Values: \n');
            %     input.prc.nativeInp
            %     fprintf('Obs Param Values: \n');
            %     input.obs.nativeInp
            %     % re-sample prc param values
            %     for j = 1:size(mod(m).prc_idx,2)
            %         input.prc.transInp(mod(m).prc_idx(j)) = ...
            %             normrnd(mod(m).prc_config.priormus(mod(m).prc_idx(j)),...
            %             abs(sqrt(mod(m).prc_config.priorsas(mod(m).prc_idx(j)))));
            %     end
            %     input.prc.nativeInp = mod(m).prc_config.transp_prc_fun(c, input.prc.transInp);
            % end
        % end

        % save sim input
        sim.sub(n,m).input = input;
        
        % Update the rng state idx
        options.rng.idx = options.rng.idx+1;
        if options.rng.idx == (length(options.rng.settings.State)+1)
            options.rng.idx = 1;
        end
    end
end

% reset rng idx
options.rng.idx = 1;


end