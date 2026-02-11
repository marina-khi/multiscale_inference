function simulate_DGP_SimBaS(t_len_, n_ts_, different_b_, pars)
%-----------------------------------------------------------------------
%  Simulate functional predictor/response pairs for SimBaS benchmarking
%-----------------------------------------------------------------------
% Inputs
%   t_len_          scalar  – # time points (T)
%   n_ts_           scalar  – # subjects  (N)
%   different_b_    row-vec – bump heights (first element 0 ⇒ size run)
%   pars            struct  – model parameters, fields listed below
%
% Required fields of pars
%   rho_            scalar  – correlation in Σbig for α
%   bump_height_    scalar  – height of null bump (size runs)
%   beta_           1×3     – regression coefficients
%   a_x_vec_        1×3     – diagonal AR(1) coeffs for VAR predictor
%   a_              scalar  – AR(1) coeff for ε
%   sigma_          scalar  – innovation s.d. for ε
%   phi_            scalar  – off-diag element in Φ for VAR innovations
%
% Output   →  CSV files simX_*.csv, simY_*.csv for every bump height
%-----------------------------------------------------------------------

% ---------- deterministic seed for full reproducibility ---------------
rng(20250805);

% ---------- covariance for random intercepts --------------------------
big_sigma = pars.rho_ * ones(n_ts_) + (1-pars.rho_) * eye(n_ts_);
alpha_vec = mvnrnd(zeros(1,n_ts_), big_sigma);   % 1×N

% ---------- pre-allocate containers -----------------------------------
error_matrix = NaN(t_len_, n_ts_);               % T×N
x_store      = NaN(t_len_, n_ts_);               % scalar exposure curves
% (y_matrices will be created per bump height)

% ---------- ARIMA model for ε -----------------------------------------
errModel = arima('AR', pars.a_, 'Constant', 0, 'Variance', pars.sigma_^2);

% ---------- VAR(1) coefficient & innovation cov -----------------------
A = diag(pars.a_x_vec_);
Phi          = pars.phi_ * ones(3);
Phi(1:4:end) = 1;    % set diagonal to 1

% ---------- helper: bump shape on [0,1] -------------------------------
bump = @(u) (u>=0 & u<=1) .* sin(pi*u);

% ---------- loop over subjects ----------------------------------------
for i = 1:n_ts_
    % --- ε_i(t) --------------------------------------------------------
    error_matrix(:,i) = simulate(errModel, t_len_);
    
    % --- latent 3-D VAR(1) predictor ----------------------------------
    nu      = mvnrnd(zeros(3,1), Phi, t_len_+10);    % (T+10)×3
    x_tmp   = zeros(3, t_len_+10);
    for tt = 2:(t_len_+10)
        x_tmp(:,tt) = A * x_tmp(:,tt-1) + nu(tt,:)';
    end
    x_mat   = x_tmp(:, 11:end)';                     % T×3
    x_store(:,i) = x_mat * pars.beta_';              % T×1   → scalar exposure
end

% ---------- iterate over bump heights ---------------------------------
grid_u = (1:t_len_)' / t_len_;                       % T×1 scaled grid
for b = different_b_
    
    % --- mean structure m_i(t) ----------------------------------------
    if b == 0         % size simulations
        m_matrix = bump(grid_u) * pars.bump_height_;          % T×1
        m_matrix = repmat(m_matrix, 1, n_ts_);                % T×N
    else               % power simulations
        m_matrix = zeros(t_len_, n_ts_);
        m_matrix(:,1) = bump(grid_u) * b;                     % only 1st curve
    end
    
    % --- build response ------------------------------------------------
    Y = alpha_vec + m_matrix + x_store + error_matrix;        % T×N
    Y = Y.';                           % transpose → N×T for SimBaS
    X = x_store.';                     % N×T exposure
    
    % --- write to disk -------------------------------------------------
    fnameX = sprintf('simX_N%d_T%d_b%g.csv', n_ts_, t_len_, b);
    fnameY = sprintf('simY_N%d_S%d_b%g.csv', n_ts_, t_len_, b);
    writematrix(X, fnameX);         %#ok<*WRIT1>
    writematrix(Y, fnameY);
    
    % --- optional: save everything in a .mat file ----------------------
    save(sprintf('DGP_b%g.mat', b), 'X', 'Y', 'alpha_vec', ...
         'm_matrix', 'error_matrix', 'pars');
end
end