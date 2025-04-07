function [vs_mod, ep_mod] = metac_est_ep(dat, ds_mod)

% mod = dat.mod;

%% load results from model inversion
ep_mod = struct();
for n = 1:size(dat.u_pe,2)
    for m = 1:size(ds_mod,2)
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n,m);
        % mod(m).sub(n).est = load(fullfile(saveDir, 'results',...
        %     'pilots', ['sub', num2str(n)], ['est_mod', num2str(m)]));
        
        % prc model
        for j = 1:size(ds_mod(m).prc_idx,2)
            ep_mod(m).prc_est(n,j) = ds_mod(m).sub(n).est.p_prc.ptrans(ds_mod(m).prc_idx(j));
            if n == 1
                ep_mod(m).prc_priormus(j) = ds_mod(m).prc_config.priormus(ds_mod(m).prc_idx(j));
                ep_mod(m).prc_priorsas(j) = ds_mod(m).prc_config.priorsas(ds_mod(m).prc_idx(j));
            end
        end
        
        %obs model
        for k = 1:size(ds_mod(m).obs_idx,2)
            ep_mod(m).obs_est(n,k) = ds_mod(m).sub(n).est.p_obs.ptrans(ds_mod(m).obs_idx(k));
            if n == 1
                ep_mod(m).obs_priormus(k) = ds_mod(m).obs_config.priormus(ds_mod(m).obs_idx(k));
                ep_mod(m).obs_priorsas(k) = ds_mod(m).obs_config.priorsas(ds_mod(m).obs_idx(k));
            end
        end           
    end
end

%% estimate pilot priors
for m = 1:size(ds_mod, 2)
    % % prc
    % model(m).prc_robmean = NaN(size(model(m).prc_priormus));
    % model(m).prc_robvar = NaN(size(model(m).prc_priorsas));
    % for j = 1:size(mod(m).prc_idx,2)
    %     [model(m).prc_robvar(1,j), model(m).prc_robmean(1,j)] = robustcov(model(m).prc_est(:,j));
    % end
    
    % obs
    ep_mod(m).obs_robmean = NaN(size(ep_mod(m).obs_priormus));
    ep_mod(m).obs_robvar = NaN(size(ep_mod(m).obs_priorsas));
    for k = 1:size(ds_mod(m).obs_idx,2)
        [ep_mod(m).obs_robvar(1,k), ep_mod(m).obs_robmean(1,k)] = robustcov(ep_mod(m).obs_est(:,k));
    end
end

%% create model space for main analysis (pilot priors)

% init model space
vs_mod = metac_create_model_space(1); %1=logit space

% set pilot priors
for m = 1:size(vs_mod, 2)
    % % prc
    % for j = 1:size(mod(m).prc_idx,2)
    %     main.mod(m).prc_config.priormus(main.mod(m).prc_idx(j)) = ...
    %         round(model(m).prc_robmean(j), 4);
    %     main.mod(m).prc_config.priorsas(main.mod(m).prc_idx(j)) = ...
    %         round(model(m).prc_robvar(j), 4);
    % end
    % main.mod(m).prc_config = tapas_align_priors_fields(...
    %     main.mod(m).prc_config);
    % obs
    for k = 1:size(ds_mod(m).obs_idx,2)
        vs_mod(m).obs_config.priormus(vs_mod(m).obs_idx(k)) = ...
            round(ep_mod(m).obs_robmean(k), 4);
        vs_mod(m).obs_config.priorsas(vs_mod(m).obs_idx(k)) = ...
            round(ep_mod(m).obs_robvar(k), 4);
    end
    vs_mod(m).obs_config = tapas_align_priors_fields(...
        vs_mod(m).obs_config);
end

end