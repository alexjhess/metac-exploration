function test_mc_pe_res_calc(dat, s)


%% get data

n = size(dat.u_bin,1);

sub = s;
x_mc = dat.u_pe(:,sub);
u_mc = dat.u_bin(:,sub);

ga = 0.65; %c.logitgamu = tapas_logit(0.65,1);
y0 = 0.5; %-0.04
w = 0.2; 



%% Giulia's implementation
% Calculate predicted response
A = zeros(n); %stores all factors in a matrix so that yhat is given by the sum of columns
for k=1:n
    for j=1:n
        if k >= j
            % A(j,k) = ga^(k-j)*x_mc(j) + 1/k*ga^k*y0 + w*(ga^(k-j)*u_mc(j) + 1/k*ga^k*y0); % GC
            A(j,k) = ga^(k-j)*x_mc(j) + 1/k*ga^k*y0 + w*(ga^(k-j)*u_mc(j)); 
        end
    end
end

z_gc = sum(A)';

%% Alex H implementation (pptx)
z_ah = zeros(n,1);
for k=1:n
    pe = zeros(n,1);
    r = zeros(n,1);
    for j = 1:n
        if k >= j
            pe(j) = ga^(k-j)*x_mc(j);
            r(j) = ga^(k-j)*u_mc(j);
            % pe(j) = ga^(k-j)*x_mc(j) + 1/k*ga^k*y0; % GC
            % r(j) = ga^(k-j)*u_mc(j) + 1/k*ga^k*y0; % GC
        end
    end
    % z_ah(k) = sum(pe) + w*sum(r); % GC
    z_ah(k) = sum(pe) + w*sum(r) + ga^k*y0;
end


%% plot comparison
figure
plot(z_gc)
hold on;
plot(z_ah, 'x')
legend('Giulia', 'Alex')

%% test results

% for k = 1:n
%     k
%     assert(z_gc(k) == z_ah(k))
% end

assert(sum(z_gc) == sum(z_ah));


end