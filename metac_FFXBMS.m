function [sumLME, pp, GBF, ABF] = metac_FFXBMS(LME)
% author: Alex J. Hess
% date: 2025-02-13
% version: 1.0
%--------------------------------------------------------------------------
% This function returns a vector containing the sum of the LME for each
% model used to generate data. Additionally, one matrix containing the
% Group Bayes Factors (GBF) and one containing the posterior model
% probabilities for every model are returned.
%
% Input:
%   - LME       n x m matrix (n=sub, m=model)
%
% Output:
%   - sumLME    ..
%   - pp        ..
%   - GBF       ..
% 
%==========================================================================


%% compute sum of the approx. LME
sumLME = sum(LME,1);
GBF =  [];

if size(LME,2) == 4
    % compute Group Bayes Factor (m4 vs m1)
    GBF = exp(sumLME(4) - sumLME(1));
    % compute average BF (ABF m4 vs m1)
    ABF = GBF^(1/size(LME,1));
    % compute the posterior model probabilities
    sumLME_maxsubtracted = sumLME - max(sumLME);
    pp = exp(sumLME_maxsubtracted)./sum(exp(sumLME_maxsubtracted));
elseif size(LME,2) == 6
    % family-level comparison
    K = 2; % nFam
    N_1 = [0 1 0 0 0 0]; % models in Family 1
    N_2 = [1 0 1 1 1 1]; % models in Family 2
    f1_priors = 1/(K*sum(N_1));
    f2_priors = 1/(K*sum(N_2));
    
    %%
    sumpp = (N_1*f1_priors+N_2*f2_priors) .* exp(sumLME);
    
    %%
    pp(1) = sum((exp(sumLME).*N_1).*(f1_priors*N_1) / sum(sumpp));
    %%
    pp(2) = sum((exp(sumLME).*N_2).*(f2_priors*N_2) / sum(sumpp));
    %%
else
    % compute the posterior model probabilities
    sumLME_maxsubtracted = sumLME - max(sumLME);
    pp = exp(sumLME_maxsubtracted)./sum(exp(sumLME_maxsubtracted));
end
   

end
