function pdat = metac_preproc(dat)


%% % Shrink all VAS responses toward 1/2 by a factor of 0.95 (s.t. in *open* unit interval)
y_pred = 0.95.*(dat.y_pred-0.5)+0.5;
y_mc = 0.95.*(dat.y_mc-0.5)+0.5;
y_c = 0.95.*(dat.y_c-0.5)+0.5;
y_tol = 0.95.*(dat.y_tol-0.5)+0.5;
y_av = 0.95.*(dat.y_av-0.5)+0.5;

%% calc PE
u_pe = dat.u_bin - y_pred;
u_pe_sq = u_pe.^2;
u_pe_pos = NaN(size(u_pe));
u_pe_pos(u_pe>=0) = u_pe(u_pe>=0);
u_pe_neg = NaN(size(u_pe));
u_pe_neg(u_pe<0) = abs(u_pe(u_pe<0));

%% create pdat struct
pdat.y_pred = y_pred;
pdat.y_mc = y_mc;
pdat.y_c = y_c;
pdat.y_tol = y_tol;
pdat.y_av = y_av;
pdat.u_pe = u_pe;
pdat.u_pe_sq = u_pe_sq;
pdat.u_pe_pos = u_pe_pos;
pdat.u_pe_neg = u_pe_neg;


end