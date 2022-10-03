clear all

% -------------------------------------------------------------------------- %
%cc%
foldername = './data2/';
miles = 75; % set it to 75, 100, or 125 miles

% -------------------------------------------------------------------------- %
% load defactored variables
load(sprintf('%sdata_main_defactored_wide_W%03d.mat', foldername, miles), 'v_id', 'v_time', 'v_reg', 'm_hp_dfac', 'm_pp_dfac', 'm_ic_dfac');
m_hp = m_hp_dfac; clear m_hp_dfac
m_pp = m_pp_dfac; clear m_pp_dfac
m_ic = m_ic_dfac; clear m_ic_dfac

% -------------------------------------------------------------------------- %
% load weight matrix W
if miles ~= 0
     m_W = csvread(sprintf('./data1/yang/W%d.csv', miles)); % (377,377)
else
     m_W = csvread('./data2/W000.csv'); % (377,377)
end

% check that every unit has at least one neighbour
v_row_sum = sum(m_W, 2); % (377,1)
N_without = sum(v_row_sum == 0);
if N_without > 0
     disp(sprintf('This weight matrix W has %d units without neighbours.', N_without))

     % locate units without neighbours
     v_ind_without = find(v_row_sum == 0);
     
     % locate units with neighbours
     v_ind_with = find(sum(m_W, 2) ~= 0);
     m_W_temp = m_W(v_ind_with, v_ind_with);
     m_W = m_W_temp;

     % remove units without neighbours
     v_id (v_ind_without) = [];
     v_reg(v_ind_without) = [];

     m_hp(v_ind_without, :) = [];
     m_pp(v_ind_without, :) = [];
     m_ic(v_ind_without, :) = [];
elseif N_without == 0
     v_ind_without = [];
end

[N, T] = size(m_hp);
disp(sprintf('   (N, T) = (%d, %d)', N, T))

% generate Wy, and Wx
m_hs = m_W * m_hp; % (N,T)
m_ps = m_W * m_pp;
m_is = m_W * m_ic;

% -------------------------------------------------------------------------- %
% build explanatory variables
K = 5; %cc%
a_x = ones(N, T - 1, K);
a_x(:, :, 2) = m_pp(:, 2:T);
a_x(:, :, 3) = m_ic(:, 2:T);
a_x(:, :, 4) = m_hs(:, 1:(T - 1)); % Wy1
a_x(:, :, 5) = m_hp(:, 1:(T - 1)); % y1

% adjust for a period lost
m_hp(:, 1) = [];
T = size(m_hp, 2);

% -------------------------------------------------------------------------- %
% define parameter space in X
%% upper bounds
m_ub_x = NaN(N, K); 
m_ub_x(:, 1) = Inf(N, 1); % intercept
m_ub_x(:, 2) = Inf(N, 1); % pp
m_ub_x(:, 3) = Inf(N, 1); % ic
m_ub_x(:, 4) = 0.995 * ones(N, 1); % Wy1
m_ub_x(:, 5) = 0.995 * ones(N, 1); % y1

%% lower bounds
m_lb_x = -m_ub_x;

% define parameter space for the psi's and sgmsq's
v_ub_Wy    = 0.995 * ones(N, 1);
v_lb_Wy    = -v_ub_Wy;
v_ub_sgmsq = Inf(N, 1);
v_lb_sgmsq = zeros(N, 1);

% put all bounds together into a single array of order (N,K+2,2)
m_lb = [v_lb_Wy m_lb_x v_lb_sgmsq];
m_ub = [v_ub_Wy m_ub_x v_ub_sgmsq];
a_b = NaN(N, K + 2, 2);
a_b(:, : , 1) = m_lb;
a_b(:, : , 2) = m_ub;

% -------------------------------------------------------------------------- %
% initial values for the optimisation procedure
v_theta_ini_tr = zeros(1, K + 2); % (1,K+2)
v_theta_ini_tr(end) = 1; % sgmsq
m_theta_ini = repmat(v_theta_ini_tr, [N 1]); % (N,K+2)

% -------------------------------------------------------------------------- %
% esitmate HSAR
results = fn_ml_Npsi_NKbeta_Nsgmsq(m_hp, a_x, m_W, a_b, m_theta_ini);

%%% check variance
%%if sum(results.m_sandwich(:) < 0) > 0
%%     error('negative variance')
%%end

% -------------------------------------------------------------------------- %
% generate residuals
[N T] = size(m_hp);
m_e = zeros(N, T);
for ii = 1:N
     psi_hat = results.m_theta(ii, 1); % (1,1)
     v_beta_hat = results.m_theta(ii, 2:(end - 1))'; % (K,1)
     for t = 1:T
          y = m_hp(ii, t);
          v_x = squeeze(a_x(ii, t, :)); % (K,1)
          y_hat = v_beta_hat' * v_x;
          m_e(ii, t) = y - y_hat;
     end
end

% -------------------------------------------------------------------------- %
% save results
save(sprintf('%sestimates_W%03d.mat', foldername, miles), ...
     'v_id', ...
     'v_reg', ...
     'results', ...
     'v_time', ... included for completeness, not needed
     'm_e', ...
     'm_W' ...
     )
