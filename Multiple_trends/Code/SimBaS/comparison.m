%======================================================================
% SimBaS_constant_covariate_FixedJitter.m   (error-free)
%======================================================================
clear; clc;

%% ----- Simulation design --------------------------------------------
N = 200;     T = 100;     V = 128;
rng(1);

t_grid    = linspace(0,1,T);
beta_true = 2*sin(2*pi*t_grid);
sigma_eps = 0.5;

X = ones(N,V) + 1e-8*randn(N,V);   % *** tiny jitter to avoid sd = 0 ***
Y = ones(N,1)*beta_true + sigma_eps*randn(N,T);

%% ----- Wavelet specs (unchanged) ------------------------------------
wavespec.comp_level = 0.99;
wavespec.alph       = [0.9 0.95 0.99 0.995 0.999 0.9999];
wavespec.nlevels    = 6;
wavespec.wavelet    = 'db4';
wavespec.wtmode     = 'zpd';
wavespec.ndim       = 1;

wavespecsy          = wavespec;
wavespecsy.xticks   = 1:T;

wavespecsx          = wavespec;
wavespecsx.xticks   = 1:V;

%% ----- Optional wPC pipeline ----------------------------------------
pcaspecsx.pca  = 1;
pcaspecsx.npcs = 1;

%% ----- Remaining inputs (as before) ---------------------------------
SIMspecs.nScalarCov = 0;     FDR = [0.01 0.05];
TSspecs.twosample   = 0;     Z   = [];

MCMCspecs = struct('B',1000,'burnin',1000,'thin',1,'propsdTheta',1,...
                   'nj_nosmooth',1,'minp',1e-14,'maxO',1e20,...
                   'minVC',1e-20,'VC0_thresh',1e-4,'delta_theta',1e-4,...
                   'thetaMLE_maxiter',1e5,'EmpBayes_maxiter',1e5,...
                   'time_update',100,'tau_prior_var',1e3,'tau_prior_idx',1,...
                   'PI_prior_var',0.06,'pi_update',1,'tau_update',1);
GBF = 0;

%% ----- Run WFFR ------------------------------------------------------
tic;
results = flmm_compress(Y, X, Z, FDR, SIMspecs, TSspecs, ...
                        MCMCspecs, wavespecsx, wavespecsy, GBF, 0, pcaspecsx);
runtime = toc;

SimBaS = results.SimBaS;
save('SimBaS_constant_covariate_FixedJitter.mat','results','SimBaS',...
     'beta_true','t_grid','runtime');

figure;
plot(t_grid, SimBaS, 'LineWidth',1.2); hold on;
yline(0.05,'--r','\alpha = 0.05','LineWidth',1);
xlabel('t'); ylabel('P_{SimBaS}(t)');
title('SimBaS – constant predictor (jitter fix)'); grid on;
