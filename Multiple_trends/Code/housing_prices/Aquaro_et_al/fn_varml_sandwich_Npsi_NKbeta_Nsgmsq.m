function [m_variance,m_sandwich,v_cov_psi0i_psi1i] = fn_varml_sandwich_Npsi_NKbeta_Nsgmsq(m_theta, m_y, m_ys, a_x, m_W)
[N T K] = size(a_x);

v_psi   = m_theta(:, 1);
m_beta  = m_theta(:, 2:(K + 1));
v_sgmsq = m_theta(:, K + 2);

% generate a_beta in order to compute the residuals
m_beta_tr = m_beta'; % K x N: (v_beta_1:...:v_beta_N)
a_beta_tr = repmat(m_beta_tr, [1 1 T]); % K x N x T
a_beta = permute(a_beta_tr, [2 3 1]); % N x T x K

v_sgm4h = v_sgmsq.^2;
v_sgm6h = v_sgmsq.^3;

m_psi  = repmat(v_psi, 1, T);
% compute residuals
m_beta_times_x = sum(a_beta .* a_x, 3); % N x T with generic element \vbeta_{i}'\vx_{it}
m_eps = m_y - (m_psi .* m_ys) - m_beta_times_x; % N x T
m_epssq = m_eps.^2;
v_ssr = sum(m_epssq, 2); % (N x 1) sum of squared residuals: sum_t (eps_it)^2

m_Psi = diag(v_psi); % N x N
m_A = eye(N) - (m_Psi * m_W);

% first derivative
m_Q = m_W * inv(m_A);

% Inverse H matrix
m_H11 = (m_Q .* m_Q') + diag(sum(m_ys.^2, 2) ./ v_sgmsq ./ T); % N x N
m_H13 = diag(sum(m_ys .* m_eps, 2) ./ v_sgm4h ./ T); % N x N
m_H33 = diag(-(1 ./ 2 ./ v_sgm4h) + (v_ssr ./ v_sgm6h ./ T)); % N x N

m_H12 = zeros(N, (N * K));
invH22 = zeros(N * K);
m_H23 = zeros((N * K), N);

for i = 1:N
     ind = ((i - 1) * K + 1):(i * K);
     v_ysi = m_ys(i, :)';
     m_Xi = permute(a_x(i, :, :), [2 3 1]); % T x K x 1 (for i-th unit)
     v_epsi = m_eps(i, :)';

     sgmsqi = v_sgmsq(i, 1);
     sgm4hi = v_sgm4h(i, 1);
     m_H12(i, ind) = (v_ysi' * m_Xi) ./ sgmsqi ./ T; % 1 x K
     invH22(ind, ind) = inv(m_Xi' * m_Xi) .* sgmsqi .* T; % K x K
     m_H23(ind, i) = (m_Xi' * v_epsi) ./ sgm4hi ./ T; % K x 1
end

m_Z11 = m_H11;
m_Z12 = [m_H12 m_H13];

invZ22 = fn_inv_partitioned_a(invH22, m_H23, m_H23', m_H33);
invH = fn_inv_partitioned_b(m_Z11, m_Z12, m_Z12', invZ22);

% J matrix
v_q = diag(m_Q);
m_sgmsq = repmat(v_sgmsq, 1, T); % N x T
m_sgm4h = repmat(v_sgm4h, 1, T); % N x T
m_dlogft_dvpsi = (m_ys .* m_eps ./ m_sgmsq) - repmat(v_q, 1, T); % N x T
v_dlogft_dvsgmsq = (m_epssq ./ m_sgm4h ./ 2) - (1 ./ 2 ./ m_sgmsq); % N x T

a_sgmsq = repmat(v_sgmsq, [1 T K]); % N x T x K
a_eps = repmat(m_eps, [1 1 K]); % N x T x K
a_dlogft_dvbeta = (a_eps .* a_x) ./ a_sgmsq; % N x T x K
a_dlogft_dvbeta_perm = permute(a_dlogft_dvbeta, [3 2 1]); % K x T x N
% collapse to an NK x T matrix as explained below:
% m_beta = [m_beta_1'; 
%           m_beta_2'; 
%           ...;
%           m_beta_N']; 
m_dlogft_dvbeta = zeros(K * N, T); % KN x T
for i = 1:N
     v_ind = ((i - 1) * K + 1):(i * K);
     m_dlogft_dvbeta(v_ind, :) = a_dlogft_dvbeta_perm(:, :, i); % K x T
end
% (N + KN + N) x T
m_dlogft_dvtheta = [m_dlogft_dvpsi; 
                    m_dlogft_dvbeta;
                    v_dlogft_dvsgmsq];
m_J = (m_dlogft_dvtheta * m_dlogft_dvtheta') ./ T;

% standard variance
v_var = diag(invH) ./ T;
m_var = zeros(N, K + 2);
m_var(:, 1) = v_var(1:N, 1);
m_var(:, 2:(K + 1)) = reshape(v_var((N + 1):(N + (K * N)), 1), K, N)';
m_var(:, K + 2) = v_var((N + (K * N) + 1):(N + (K * N) + N), 1);
m_variance = m_var;

% sandwich variance
m_invH_J_invH = invH * m_J * invH;
v_var = diag(m_invH_J_invH) ./ T;
m_var = zeros(N, K + 2);
m_var(:, 1) = v_var(1:N, 1);
m_var(:, 2:(K + 1)) = reshape(v_var((N + 1):(N + (K * N)), 1), K, N)';
m_var(:, K + 2) = v_var((N + (K * N) + 1):(N + (K * N) + N), 1);
m_sandwich = m_var;

% sandwich cov(psi0i, psi1i) %cc% !! THIS BIT OF CODE IS VALID ONLY WITH MODEL 4 (Wy0+x, where x=1+pp+ic+Wy1+y1)
m_varcov = m_invH_J_invH / T;
m_varcov_12 = m_varcov(1:N, (N + 1):(N + (K * N))); % (N,NK)
v_cov_psi0i_psi1i = NaN(N, 1);
for ii = 1:N
     % take entries (1,4), (2,9), (3,14), (4,19) corresponding to cov(psi0,1, psi1,1), cov(psi0,2, psi1,2), cov(psi0,3, psi1,3), etc.
     v_cov_psi0i_psi1i(ii, 1) = m_varcov_12(ii, (ii * 5) - 1); %cc%
end
end

% ----------------------------------------------------------------------- %
function invH = fn_inv_partitioned_a(invA, m_B, m_C, m_D)
m_C_invA = m_C * invA;
m_E = m_D - (m_C_invA * m_B);
invE = inv(m_E);
m_invA_B_invE = invA * m_B * invE;

invH11 = invA + (m_invA_B_invE * m_C_invA);
invH12 = -m_invA_B_invE;
invH21 = -invE * m_C_invA;
invH22 = invE;

invH = [[invH11 invH12]; ...
        [invH21 invH22]];
end

% ----------------------------------------------------------------------- %
function invH = fn_inv_partitioned_b(m_A, m_B, m_C, invD)
m_B_invD = m_B * invD;
m_F = m_A - (m_B_invD * m_C);
invF = inv(m_F);
m_invD_C_invF = invD * m_C * invF;

invH11 = invF;
invH12 = -invF * m_B_invD;
invH21 = -m_invD_C_invF;
invH22 = invD + (m_invD_C_invF * m_B_invD);

invH = [[invH11 invH12]; ...
        [invH21 invH22]];
end
