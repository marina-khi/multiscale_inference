function l_out = fn_ml_Npsi_NKbeta_Nsgmsq(m_y, a_x, m_W, a_b, m_theta_ini)
% ----------------------------------------------------------------------- %
% PURPOSE 
% To compute the coefficient estimates of the spatial autoregressive panel data
% models with heterogeneous coefficients (HSAR), see reference below.
%
% ----------------------------------------------------------------------- %
% REQUIRED INPUT PARAMETERS
% - m_y = NxT matrix of outcomes;
% - a_x = NxTxK array of exogenous regressors (including a column of 1s for the
%   intercept);
% - m_W = NxN spatial weight matrix;
% - m_theta_ini = Nx(K+2) matrix containing the starting values
%   for numerical optimization. This vector is organised as follows: 
%  
%   [psi; beta; sigmasq]
%  
%   where psi and sigmasq are column vectors of length N, and beta is a matrix of order NxK,
%   representing:
%     * psi: Nx1 vector of heterogeneous spatial autoregressive parameters;
%     * beta: NxK matrix of heterogeneous parameters corresping to the exogeous variables;
%     * sigmasq: Nx1 vector of heterogeous variance of the ideosincratic error;
%  
%   If m_theta_ini is not supplied, the matrix:
%  
%   [zeros(N,K+1); ones(N,1)]
%  
%   is used instead.
%
% Notes: 
% - the panel is assumed to be balanced, i.e. there are no missing observations;
% - all units have at least one neighbour, i.e. "sum(m_W, 2) ~= 0";
% - the array a_x includes an NxT matrix of ones (which corresponds to an heterogeneous intercept);
% - it is assumed that m_W is row-standarised. Different types of standardisation
%   can be easily accommodated by changing the values pm0.995 in Matlab's fmincon functin.
%
% OUTPUT PARAMETERS
% l_out = structure containing three Nx(K+2) matrices plus a vector of length N:
% - l_out.m_theta: coefficient estimates, see note below;
% - l_out.m_variance: estimated variance of coefficient estimates using the standard formula;
% - l_out.m_sandwich: estimated variance of coefficient estimates using the sandwich formula;
% - l_out.v_cov_psi0i_psi1i: estimated variance-covariance matrix cov(psi0i,psi1i);
%
% Notes:
% - The columns of the three matrix above are organised as follows:
%
%   [psi beta sigmasq]
%
%   where beta is an NxK matrix.
%
% ----------------------------------------------------------------------- %
% REFERENCES
% Aquaro, M., Bailey, N. and Pesaran, H. P. (2020). A Quasi Maximum Likelihood
% Estimation of Spatial Models with Heterogeneous Parameters, Journal of
% Applied Econometrics, forthcoming
%
% ----------------------------------------------------------------------- %
% DISCLAIMER
%
% This code is offered with no guarantees. Please let us know if you find any
% bugs or encounter any problems while using this code. All feedback is
% appreciated.
%
% Please address any question to:
% Michele Aquaro
% Joint Research Centre
% European Commission
% Email: Michele.AQUARO at ec.europa.eu
%
% Last update: Jun 2020.
% ----------------------------------------------------------------------- %
 
[N T K] = size(a_x);
m_ys = m_W * m_y;

options = optimset('Display', 'iter-detailed', ...
                   'GradObj', 'on', ... 
                   'DerivativeCheck', 'off', ...
                   'Diagnostics', 'on', ...
                   'Hessian', 'on');
options.Algorithm = 'trust-region-reflective';

if nargin < 5
     v_theta_ini = [zeros((K + 1) * N, 1); ones(N, 1)];
else
     v_psi_ini = m_theta_ini(:, 1); % (N,1)
     v_sgmsq_ini = m_theta_ini(:, end); % (N,1) 
     m_beta_ini = m_theta_ini(:, 2:(end - 1)); % (N,K)
     m_beta_ini_tr = m_beta_ini'; % (K,N)
     v_beta_ini = m_beta_ini_tr(:); % (KN,1)
     v_theta_ini = [v_psi_ini; v_beta_ini; v_sgmsq_ini]; % ((K+2)N,1)
end

% lower/upper bounds in the constrained optimisation
%% in matrix form
m_lb = a_b(:, :, 1); % (N,K+2)
m_ub = a_b(:, :, 2);

%% psi
v_lb_psi = m_lb(:, 1); % (N,1)
v_ub_psi = m_ub(:, 1);

%% sgmsq
v_lb_sgmsq = m_lb(:, end); % (N,1)
v_ub_sgmsq = m_ub(:, end);

%% beta
m_lb_beta = m_lb(:, 2:(end - 1)); % (N,K)
m_ub_beta = m_ub(:, 2:(end - 1));
m_lb_beta_tr = m_lb_beta'; % (K,N)
m_ub_beta_tr = m_ub_beta';
v_lb_beta = m_lb_beta_tr(:); % (KN,1)
v_ub_beta = m_ub_beta_tr(:);

v_lb = [v_lb_psi; v_lb_beta; v_lb_sgmsq]; % (NP,1)
v_ub = [v_ub_psi; v_ub_beta; v_ub_sgmsq];

fn_anonym = @(v_theta)fn_uncloglk_Npsi_NKbeta_Nsgmsq(v_theta, m_y, m_ys, a_x, m_W);
[v_theta fval exitflag output lambda grad hessian] = fmincon(fn_anonym, v_theta_ini, [], [], [], [], v_lb, v_ub, [], options); exitflag
v_psi   = v_theta(1:N, 1);
v_beta  = v_theta((N + 1):(N + (K * N)), 1);
v_sgmsq = v_theta((N + (K * N) + 1):(N + (K * N) + N), 1);
m_beta_tr = reshape(v_beta, K, N);
m_theta = [v_psi m_beta_tr' v_sgmsq];

% \hat{var}(\hat{v_theta}) to be computed separately to exploit blocks when taking the inverse
[m_variance,m_sandwich,v_cov_psi0i_psi1i] = fn_varml_sandwich_Npsi_NKbeta_Nsgmsq(m_theta, m_y, m_ys, a_x, m_W);

l_out.m_theta = m_theta;
l_out.m_variance = m_variance;
l_out.m_sandwich = m_sandwich;
l_out.v_cov_psi0i_psi1i = v_cov_psi0i_psi1i;
end

% ----------------------------------------------------------------------- %
% phi(vtheta):=-(logL(vtheta) / T)
function [ofun v_gradient m_H] = fn_uncloglk_Npsi_NKbeta_Nsgmsq(v_theta, m_y, m_ys, a_x, m_W)
[N T K] = size(a_x);

v_psi   = v_theta(1:N, 1);
v_beta  = v_theta((N + 1):(N + (K * N)), 1); % KN x 1
v_sgmsq = v_theta((N + (K * N) + 1):(N + (N * K) + N), 1); % N x 1

% generate a_beta in order to compute the residuals
m_beta_tr = reshape(v_beta, K, N); % K x N: (v_beta_1:...:v_beta_N)
a_beta_tr = repmat(m_beta_tr, [1 1 T]); % K x N x T
a_beta = permute(a_beta_tr, [2 3 1]); % N x T x K

v_sgm4h = v_sgmsq.^2;
v_sgm6h = v_sgmsq.^3;

m_psi  = repmat(v_psi, 1, T);
% compute residuals
m_beta_times_x = sum(a_beta .* a_x, 3); % N x T with generic element \vbeta_{i}'\vx_{it}
m_eps = m_y - (m_psi .* m_ys) - m_beta_times_x; % N x T
v_ssr = sum(m_eps.^2, 2); % (N x 1) sum of squared residuals: sum_t (eps_it)^2
% standardize
sssr = sum(v_ssr ./ v_sgmsq, 1); % sum of squared standardized-residuals: sum_i sum_t (eps_it / sgm_i)^2

m_Psi = diag(v_psi); % N x N
m_A = eye(N) - (m_Psi * m_W);

det_mA = det(m_A);
if det_mA <= 0
     error('Michele: error in fn_uncloglk_Npsi_NKbeta_Nsgmsq')
end
constant = log(2 * pi) * N / 2;
first_part = -log(det_mA);
secon_part = sum(log(v_sgmsq)) / 2;
third_part = sssr / T / 2;
ofun = constant + first_part + secon_part + third_part;

% first derivative
m_Q = m_W * inv(m_A);
v_dphi_dvpsi = diag(m_Q) - (sum(m_ys .* m_eps, 2) ./ v_sgmsq ./ T);

m_sgmsq = repmat(v_sgmsq, [1 K]); % N x K
a_eps = repmat(m_eps, [1, 1, K]); % N x T x K
a_X_times_eps = sum(a_x .* a_eps, 2); % N x 1 x K
m_X_times_eps = squeeze(a_X_times_eps); % N x K
m_X_times_eps_divided_sgmsq = m_X_times_eps ./ m_sgmsq; % N x K
m_X_times_eps_divided_sgmsq_tr = m_X_times_eps_divided_sgmsq'; % K x N
v_dphi_dvbeta = -m_X_times_eps_divided_sgmsq_tr(:) ./ T; % NK x 1

v_dphi_dvsgmsq = (1 ./ 2 ./ v_sgmsq) - (v_ssr ./ v_sgm4h ./ T ./ 2);
v_gradient = [v_dphi_dvpsi;
              v_dphi_dvbeta;
              v_dphi_dvsgmsq];

% second derivative
m_H11 = (m_Q .* m_Q') + diag(sum(m_ys.^2, 2) ./ v_sgmsq ./ T); % N x N
m_H13 = diag(sum(m_ys .* m_eps, 2) ./ v_sgm4h ./ T); % N x N
m_H33 = diag(-(1 ./ 2 ./ v_sgm4h) + (v_ssr ./ v_sgm6h ./ T)); % N x N

m_H12 = zeros(N, (N * K));
m_H22 = zeros(N * K);
m_H23 = zeros((N * K), N);

% Note 1: the loop below can be probably eliminated in the same way v_dphi_dvbeta is computed above
% Note 2: the code can be made probably faster by using sparse matrices
for i = 1:N
     ind = ((i - 1) * K + 1):(i * K);
     v_ysi = m_ys(i, :)';
     m_Xi = permute(a_x(i, :, :), [2 3 1]); % T x K x 1 (for i-th unit)
     v_epsi = m_eps(i, :)';

     sgmsqi = v_sgmsq(i, 1);
     sgm4hi = v_sgm4h(i, 1);
     m_H12(i, ind) = (v_ysi' * m_Xi) ./ sgmsqi ./ T; % 1 x K
     m_H22(ind, ind) = (m_Xi' * m_Xi) ./ sgmsqi ./ T; % K x K
     m_H23(ind, i) = (m_Xi' * v_epsi) ./ sgm4hi ./ T; % K x 1
end

m_H = [[m_H11  m_H12  m_H13]; ...
       [m_H12' m_H22  m_H23]; ...
       [m_H13' m_H23' m_H33]];

end
