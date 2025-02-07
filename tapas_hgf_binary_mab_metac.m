function [traj, infStates] = tapas_hgf_binary_mab_metac(r, p, varargin)
% Calculates the trajectories of the agent's representations under the HGF in a multi-armed bandit
% situation with binary outcomes
%
% This function can be called in two ways:
% 
% (1) tapas_hgf_binary_mab(r, p)
%   
%     where r is the structure generated by tapas_fitModel and p is the parameter vector in native space;
%
% (2) tapas_hgf_binary_mab(r, ptrans, 'trans')
% 
%     where r is the structure generated by tapas_fitModel, ptrans is the parameter vector in
%     transformed space, and 'trans' is a flag indicating this.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2017 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans')
    p = tapas_hgf_binary_mab_transp(r, p);
end

% Number of levels
try
    l = r.c_prc.n_levels;
catch
    l = (length(p)+1)/5;
    
    if l ~= floor(l)
        error('tapas:hgf:UndetNumLevels', 'Cannot determine number of levels');
    end
end

% Number of bandits
try
    b = r.c_prc.n_bandits;
catch
    error('tapas:hgf:NumOfBanditsConfig', 'Number of bandits has to be configured in r.c_prc.n_bandits.');
end

% Coupled updating
% This is only allowed if there are 2 bandits. We here assume that the mu1hat for the two bandits
% add to unity.
coupled = false;
if r.c_prc.coupled == true
    if b == 2
        coupled = true;
    else
        error('tapas:hgf:HgfBinaryMab:CoupledOnlyForTwo', 'Coupled updating can only be configured for 2 bandits.');
    end
end

% Unpack parameters
mu_0 = p(1:l);
sa_0 = p(l+1:2*l);
rho  = p(2*l+1:3*l);
ka   = p(3*l+1:4*l-1);
om   = p(4*l:5*l-2);
th   = exp(p(5*l-1));

% Add dummy "zeroth" trial
u = [0; r.u(:,1)];
try % For estimation
%     u_mab = [1; r.y(:,1)];
%     r.irr = r.irr;
% catch % For simulation
    u_mab = [1; r.u(:,2)];
    r.irr = find(isnan(r.u(:,2)));
end


% Number of trials (including prior)
n = size(u,1);

% Construct time axis
if r.c_prc.irregular_intervals
    if size(u,2) > 1
        t = [0; r.u(:,end)];
    else
        error('tapas:hgf:InputSingleColumn', 'Input matrix must contain more than one column if irregular_intervals is set to true.');
    end
else
    t = ones(n,1);
end

% Initialize updated quantities

% Representations
mu = NaN(n,l,b);
pi = NaN(n,l,b);

% Other quantities
muhat = NaN(n,l,b);
pihat = NaN(n,l,b);
v     = NaN(n,l);
w     = NaN(n,l-1);
da    = NaN(n,l);

% Representation priors
% Note: first entries of the other quantities remain
% NaN because they are undefined and are thrown away
% at the end; their presence simply leads to consistent
% trial indices.
mu(1,1,:) = tapas_sgm(mu_0(2), 1);
pi(1,1,:) = Inf;
mu(1,2:end,:) = repmat(mu_0(2:end),[1 1 b]);
pi(1,2:end,:) = repmat(1./sa_0(2:end),[1 1 b]);

% Pass through representation update loop
for k = 2:1:n
    if not(ismember(k-1, r.ign))
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
        % 2nd level prediction
        muhat(k,2,:) = mu(k-1,2,:) +t(k) *rho(2);
        
        % 1st level
        % ~~~~~~~~~
        % Prediction
        muhat(k,1,:) = tapas_sgm(ka(1) *muhat(k,2,:), 1);
        
        % Precision of prediction
        pihat(k,1,:) = 1/(muhat(k,1,:).*(1 -muhat(k,1,:)));

        % Updates
        pi(k,1,:) = pihat(k,1,:);
        pi(k,1,u_mab(k)) = Inf;
        
        mu(k,1,:) = muhat(k,1,:);
        mu(k,1,u_mab(k)) = u(k);

        % Prediction error
        da(k,1) = mu(k,1,u_mab(k)) -muhat(k,1,u_mab(k));

        % 2nd level
        % ~~~~~~~~~
        % Prediction: see above
        
        % Precision of prediction
        pihat(k,2,:) = 1/(1/pi(k-1,2,:) +exp(ka(2) *mu(k-1,3,:) +om(2)));

        % Updates
        pi(k,2,:) = pihat(k,2,:);
        pi(k,2,u_mab(k)) = pihat(k,2,u_mab(k)) +ka(1)^2/pihat(k,1,u_mab(k));

        mu(k,2,:) = muhat(k,2,:);
        mu(k,2,u_mab(k)) = muhat(k,2,u_mab(k)) +ka(1)/pi(k,2,u_mab(k)) *da(k,1);

        % Volatility prediction error
        da(k,2) = (1/pi(k,2,u_mab(k)) +(mu(k,2,u_mab(k)) -muhat(k,2,u_mab(k)))^2) *pihat(k,2,u_mab(k)) -1;

        if l > 3
            % Pass through higher levels
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~
            for j = 3:l-1
                % Prediction
                muhat(k,j,:) = mu(k-1,j,:) +t(k) *rho(j);
                
                % Precision of prediction
                pihat(k,j,:) = 1/(1/pi(k-1,j,:) +t(k) *exp(ka(j) *mu(k-1,j+1,:) +om(j)));

                % Weighting factor
                v(k,j-1) = t(k) *exp(ka(j-1) *mu(k-1,j,u_mab(k)) +om(j-1));
                w(k,j-1) = v(k,j-1) *pihat(k,j-1,u_mab(k));

                % Updates
                pi(k,j,:) = pihat(k,j,:) +1/2 *ka(j-1)^2 *w(k,j-1) *(w(k,j-1) +(2 *w(k,j-1) -1) *da(k,j-1));

                if pi(k,j,1) <= 0
                    error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                end

                mu(k,j,:) = muhat(k,j,:) +1/2 *1/pi(k,j) *ka(j-1) *w(k,j-1) *da(k,j-1);
    
                % Volatility prediction error
                da(k,j) = (1/pi(k,j,u_mab(k)) +(mu(k,j,u_mab(k)) -muhat(k,j,u_mab(k)))^2) *pihat(k,j,u_mab(k)) -1;
            end
        end

        % Last level
        % ~~~~~~~~~~
        % Prediction
        muhat(k,l,:) = mu(k-1,l,:) +t(k) *rho(l);
        
        % Precision of prediction
        pihat(k,l,:) = 1/(1/pi(k-1,l,:) +t(k) *th);

        % Weighting factor
        v(k,l)   = t(k) *th;
        v(k,l-1) = t(k) *exp(ka(l-1) *mu(k-1,l,u_mab(k)) +om(l-1));
        w(k,l-1) = v(k,l-1) *pihat(k,l-1,u_mab(k));
        
        % Updates
        pi(k,l,:) = pihat(k,l,:) +1/2 *ka(l-1)^2 *w(k,l-1) *(w(k,l-1) +(2 *w(k,l-1) -1) *da(k,l-1));
 
        if pi(k,l,1) <= 0
            error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
        end

        mu(k,l,:) = muhat(k,l,:) +1/2 *1/pi(k,l,:) *ka(l-1) *w(k,l-1) *da(k,l-1);
    
        % Volatility prediction error
        da(k,l) = (1/pi(k,l,u_mab(k)) +(mu(k,l,u_mab(k)) -muhat(k,l,u_mab(k)))^2) *pihat(k,l,u_mab(k)) -1;
        
        if coupled == true
            if u_mab(k) == 1
                mu(k,1,2) = 1 -mu(k,1,1);
                mu(k,2,2) = tapas_logit(1 -tapas_sgm(mu(k,2,1), 1), 1);
            elseif u_mab(k) == 2
                mu(k,1,1) = 1 -mu(k,1,2);
                mu(k,2,1) = tapas_logit(1 -tapas_sgm(mu(k,2,2), 1), 1);
            end
        end
    else

        mu(k,:,:) = mu(k-1,:,:); 
        pi(k,:,:) = pi(k-1,:,:);

        muhat(k,:,:) = muhat(k-1,:,:);
        pihat(k,:,:) = pihat(k-1,:,:);
        
        v(k,:)  = v(k-1,:);
        w(k,:)  = w(k-1,:);
        da(k,:) = da(k-1,:);
        
    end
end

% Remove representation priors
mu(1,:,:)  = [];
pi(1,:,:)  = [];

% Check validity of trajectories
if any(isnan(mu(:))) || any(isnan(pi(:)))
    error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
else
    % Check for implausible jumps in trajectories
    dmu = diff(mu(:,2:end));
    dpi = diff(pi(:,2:end));
    rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
    rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

    jumpTol = 16;
    if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
    end
end

% Remove other dummy initial values
muhat(1,:,:) = [];
pihat(1,:,:) = [];
v(1,:)       = [];
w(1,:)       = [];
da(1,:)      = [];
u_mab(1)         = [];

% Responses on regular trials
yreg = u_mab;
yreg(r.irr) =[];

% Implied learning rate at the first level
mu2          = squeeze(mu(:,2,:));
mu2(r.irr,:) = [];
mu2obs       = mu2(sub2ind(size(mu2), (1:size(mu2,1))', yreg));

mu1hat          = squeeze(muhat(:,1,:));
mu1hat(r.irr,:) = [];
mu1hatobs       = mu1hat(sub2ind(size(mu1hat), (1:size(mu1hat,1))', yreg));

upd1 = tapas_sgm(ka(1)*mu2obs,1) -mu1hatobs;

dareg          = da;
dareg(r.irr,:) = [];

lr1reg = upd1./dareg(:,1);
lr1    = NaN(n-1,1);
lr1(setdiff(1:n-1, r.irr)) = lr1reg;

% Create result data structure
traj = struct;

traj.mu     = mu;
traj.sa     = 1./pi;

traj.muhat  = muhat;
traj.sahat  = 1./pihat;

traj.v      = v;
traj.w      = w;
traj.da     = da;

% Updates with respect to prediction
traj.ud = mu -muhat;

% Psi (precision weights on prediction errors)
psi = NaN(n-1,l);

pi2          = squeeze(pi(:,2,:));
pi2(r.irr,:) = [];
pi2obs       = pi2(sub2ind(size(pi2), (1:size(pi2,1))', yreg));

psi(setdiff(1:n-1, r.irr), 2) = 1./pi2obs;

for i=3:l
    pihati          = squeeze(pihat(:,i-1,:));
    pihati(r.irr,:) = [];
    pihatiobs       = pihati(sub2ind(size(pihati), (1:size(pihati,1))', yreg));
    
    pii          = squeeze(pi(:,i,:));
    pii(r.irr,:) = [];
    piiobs       = pii(sub2ind(size(pii), (1:size(pii,1))', yreg));
    
    psi(setdiff(1:n-1, r.irr), i) = pihatiobs./piiobs;
end

traj.psi = psi;

% Epsilons (precision-weighted prediction errors)
epsi        = NaN(n-1,l);
epsi(:,2:l) = psi(:,2:l) .*da(:,1:l-1);
traj.epsi   = epsi;

% Full learning rate (full weights on prediction errors)
wt        = NaN(n-1,l);
wt(:,1)   = lr1;
wt(:,2)   = psi(:,2);
wt(:,3:l) = 1/2 *(v(:,2:l-1) *diag(ka(2:l-1))) .*psi(:,3:l);
traj.wt   = wt;

% Create matrices for use by the observation model
infStates = NaN(n-1,l,b,4);
infStates(:,:,:,1) = traj.muhat;
infStates(:,:,:,2) = traj.sahat;
infStates(:,:,:,3) = traj.mu;
infStates(:,:,:,4) = traj.sa;

return;
