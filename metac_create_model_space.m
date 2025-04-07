function [mod, bo] = metac_create_model_space(mc_only)
% [ModelSpace] = metac_create_model_space()
%
% Creates and outputs a struct of the model space containing config files 
% for perceptual and response models.
%
% INPUT
%   argin            type           
%
%   OPTIONAL:
%
% OUTPUT    
%   ModelSpace       struct      Struct containing prc & obs model configs
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

%% init struct
mod = struct();

%% Bayes Optimal (dummy obs model)
bo = struct();

%% init model space

if nargin > 0
    
    %% Response Models
    %----------------------------------------------------------------------
    % only y_mc
    if mc_only == 1 % logit space

        % M1: Null model, logit(y)
        mod(1).obs_name = 'null logit';
        mod(1).obs_config = logit_mc_null_obs_config;
        mod(1).obs = 'logit_mc_null_obs';
    
        % M2: AR(k) model (Res), logit(y)
        mod(2).obs_name = 'AR(k) Res recoded logit';
        mod(2).obs_config = logit_mc_autoreg_res_obs_config;
        mod(2).obs = 'logit_mc_autoreg_res_obs';

        % M3: AR(k) model (PE), logit(y)
        mod(3).obs_name = 'AR(k) PEsq logit';
        mod(3).obs_config = logit_mc_autoreg_pe_obs_config;
        mod(3).obs = 'logit_mc_autoreg_pe_obs';
        
        % M4: AR(k) model (PE + RES), logit(y)
        mod(4).obs_name = 'AR(k) sqPE + Res recoded logit';
        mod(4).obs_config = logit_mc_autoreg_pe_res_obs_config;
        mod(4).obs = 'logit_mc_autoreg_pe_res_obs';

        % M5: AR(k) model (PE + RES), logit(y)
        mod(5).obs_name = 'AR(k) sqPE + Res recoded logit';
        mod(5).obs_config = logit_mc_autoreg_pe_res_obs_config;
        mod(5).obs = 'logit_mc_autoreg_pe_res_obs';

    elseif mc_only == 2 % native space
        
        % M1: Null model
        mod(1).obs_name = 'null';
        mod(1).obs_config = mc_null_obs_config;
    
        % M2: AR(k) model (Res)
        mod(2).obs_name = 'AR(k) Res recoded';
        mod(2).obs_config = mc_autoreg_res_obs_config;

        % M3: AR(k) model (PE)
        mod(3).obs_name = 'AR(k) PEsq';
        mod(3).obs_config = mc_autoreg_pe_obs_config;

        % M4: AR(k) model (PE + RES, normal)
        mod(4).obs_name = 'AR(k) PEsq + Res recoded';
        mod(4).obs_config = mc_autoreg_pe_res_obs_config;
        
        % M5: AR(k) model (PE + RES, beta)
        mod(5).obs_name = 'AR(k) beta PEsq + Res recoded';
        mod(5).obs_config = beta_mc_autoreg_pe_res_obs_config;

    end
        
    %% Perceptual Model

    if mc_only == 1 % logit space

        % raw PE squared
        for m = 1:size(mod,2)
            if m == 5
                mod(m).prc_name = 'logit_raw PE squared';
                mod(m).prc_config = logit_raw_PEsq_config;
                mod(m).prc = 'logit_raw_PEsq';
            elseif m < 5
                mod(m).prc_name = 'logit raw PE signed';
                mod(m).prc_config = logit_raw_PE_config;
                mod(m).prc = 'logit_raw_PE';
            end
        end

    elseif mc_only == 2 % native space

        % raw PE squared
        for m = 1:size(mod,2)
            mod(m).prc_name = 'raw PE squared';
            mod(m).prc_config = raw_PEsq_config;
            mod(m).prc = 'raw_PEsq';
        end

    end
    

else % multimodal obs (y_int + y_mc) 

    %% Bayes Optimal (dummy obs model)
    bo.obs_config = tapas_bayes_optimal_binary_config;
    bo.prc_config = tapas_hgf_binary_mab_2mu0_metac_config;
    bo.prc_config.n_bandits = 2;

    % find free pars
    prc_idx = bo.prc_config.priorsas;
    prc_idx(isnan(prc_idx)) = 0;
    bo.prc_idx = find(prc_idx);
    

    %% Response Models
    %--------------------------------------------------------------------------
    % combining y_int and y_mc
    
    % M1: Null model
    mod(1).obs_name = 'null';
    mod(1).obs_config = mc_null_inobs_mab_config;
    
    % M2: AR(k) model (Res)
    mod(2).obs_name = 'AR(k) Res';
    mod(2).obs_config = mc_autoreg_res_inobs_config;
    
    % M3: AR(k) model (PE)
    mod(3).obs_name = 'AR(k) PE';
    mod(3).obs_config = mc_autoreg_pe_inobs_mab_config;
    
    % M4: AR(k) model (PE + RES)
    mod(4).obs_name = 'AR(k) PE + Res (bug)';
    mod(4).obs_config = mc_autoreg_pe_res_inobs_mab_config; % bug in y_mc calc! (y0 influence)
    
    % M5: AR(k) model (PE + RES)
    mod(5).obs_name = 'AR(k) PE + Res';
    mod(5).obs_config = mc_autoreg_pe_res_mab_obs_config;
    
    % M6: AR(k) model (PE + RES)
    mod(6).obs_name = 'AR(k) sqPE + Res';
    mod(6).obs_config = mc_autoreg_sqpe_res_mab_obs_config;
    
    % M7: AR(k) model (PE + RES)
    mod(7).obs_name = 'AR(k) sqPE + Res recoded';
    mod(7).obs_config = mc_autoreg_sqpe_negres_mab_obs_config;
    
    % M8: AR(k) model (PE + RES)
    mod(8).obs_name = 'AR(k) beta sqPE + Res recoded';
    mod(8).obs_config = beta_mc_autoreg_sqpe_negres_mab_obs_config;
    
    % cannot use same model comparison, since y_mc is transformed!!
    % M9: AR(k) model (PE + RES)
    mod(9).obs_name = 'AR(k) logit sqPE + Res recoded';
    mod(9).obs_config = logit_mc_autoreg_sqpe_negres_mab_obs_config;
    
    
    %% Perceptual Model
    % 3-level binary HGF mab 2mu0
    for m = 1:size(mod,2)
        mod(m).prc_config = tapas_hgf_binary_mab_2mu0_metac_config;
        mod(m).prc_config.n_bandits = 2;
    end

end



%% find free parameters

% obs
for m = 1:size(mod,2)
    % prc
    prc_idx = mod(m).prc_config.priorsas;
    prc_idx(isnan(prc_idx)) = 0;
    mod(m).prc_idx = find(prc_idx);
    % obs
    obs_idx = mod(m).obs_config.priorsas;
    obs_idx(isnan(obs_idx)) = 0;
    mod(m).obs_idx = find(obs_idx);
end


end