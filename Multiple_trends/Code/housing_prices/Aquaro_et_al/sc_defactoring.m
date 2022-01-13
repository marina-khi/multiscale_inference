clear all

% -------------------------------------------------------------------------- %
%cc%
foldername = './data2/';
miles = 75; % set it to 75, 100, or 125 miles, or 0 for the 75-miles inverse distance as in Table~xx in the supplement
CD = true; % save results as csv for CD test in GAUSS?

% -------------------------------------------------------------------------- %
% load data
tab = readtable('./data1/yang/data_main.csv', ...
     'Delimiter', ',', ...
     'ReadVariableNames', true, ...
     'Format', '%{yyyyQQ}D %u %u %f %f %f %f'); % long format (NT,7)

v_time = tab.quarter;
v_id   = tab.msacode;
v_reg  = tab.regioncode;
v_cpi  = tab.cpi;
v_hp   = tab.hp;
v_pop  = tab.pop;
v_pipc = tab.pipc;

NT = length(v_id); % 377 * 160
N = length(unique(v_id)); % 377
T = length(unique(v_time)); % 160
if NT ~= (N * T)
     error('The panel is unbalanced')
end
disp(sprintf('(N, T) = (%d, %d)', N, T))

% recode regions according to BEA classification
v_reg_temp = zeros(NT, 1);
v_reg_temp(v_reg == 4) = 1;
v_reg_temp(v_reg == 3) = 2;
v_reg_temp(v_reg == 2) = 3;
v_reg_temp(v_reg == 5) = 4;
v_reg_temp(v_reg == 7) = 5;
v_reg_temp(v_reg == 8) = 6;
v_reg_temp(v_reg == 6) = 7;
v_reg_temp(v_reg == 1) = 8;
v_reg = v_reg_temp; clear v_reg_temp

% from long to wide
m_id   = reshape(v_id  , [N, T]); % (N,T)
m_time = reshape(v_time, [N, T]);
m_reg  = reshape(v_reg , [N, T]); 
m_cpi  = reshape(v_cpi , [N, T]);
m_hp   = reshape(v_hp  , [N, T]);
m_pop  = reshape(v_pop , [N, T]);
m_pipc = reshape(v_pipc, [N, T]);

% MSA and regional identifiers
v_id   = m_id (:, 1); % (N,1)
v_reg  = m_reg(:, 1); % (N,1) 

% clear
clear m_id m_reg
clear v_time v_cpi v_hp v_pop v_pipc

%% construct variables
% log of real house price
m_hpr_ln = reallog(m_hp ./ m_cpi);

% canges in log of real house price
m_hp = (m_hpr_ln(:, 2:T) - m_hpr_ln(:, 1:(T - 1))) * 100;

% canges in log of population
m_pop_ln = reallog(m_pop);
m_pp = (m_pop_ln(:, 2:T) - m_pop_ln(:, 1:(T - 1))) * 100;

% canges in log of real income per capita
m_pipc_ln = reallog(m_pipc ./ m_cpi);
m_ic = (m_pipc_ln(:, 2:T) - m_pipc_ln(:, 1:(T - 1))) * 100;

% adjust for loosing one time-period because of differencing
v_time = m_time(1, 2:T)'; T = length(v_time); clear m_time

% spatially lagged changes in log of real house price
if miles == 0
     m_W = csvread(sprintf('./data2/%s', 'W000.csv')); %cc%
else 
     m_W = csvread(sprintf('./data1/yang/W%d.csv', miles)); %cc%
end
m_hs = m_W * m_hp;

if miles ~= 0 && CD
     % save m_hp (before defactoring) for CD- and alpha-test (currently available only in GAUSS)
     writematrix(m_hp, sprintf('%shp_raw.txt', foldername), 'Delimiter', 'tab');
end

%%% de-factoring and de-seasonlising
%% quarterly time dummies
v_quarter = quarter(v_time); % (T,1)
v_dummy_q2 = (v_quarter == 2); % (T,1)
v_dummy_q3 = (v_quarter == 3);
v_dummy_q4 = (v_quarter == 4);
m_quarterly_dummies = [v_dummy_q2 v_dummy_q3 v_dummy_q4]; % (T,3)

%% factors
% initialise
m_f_hp = ones(T, 2); % (T,2)
m_f_pp = ones(T, 2);
m_f_ic = ones(T, 2);
m_f_hs = ones(T, 2);

%% US factors
v_f_hp_usa = mean(m_hp, 1)'; % (T,1)
v_f_pp_usa = mean(m_pp, 1)';
v_f_ic_usa = mean(m_ic, 1)';
v_f_hs_usa = mean(m_hs, 1)';

m_f_hp(:, 1) = v_f_hp_usa;
m_f_pp(:, 1) = v_f_pp_usa;
m_f_ic(:, 1) = v_f_ic_usa;
m_f_hs(:, 1) = v_f_hs_usa;

%% regional factors
% list of regions
v_reg_unique = unique(v_reg); % (R,1)
R = length(v_reg_unique);

% the code is written with regions being ordered from 1:8 - stop if not
if ~isequal(v_reg_unique, [1:8]')
     error('the code is written with regions strictly being ordered from 1:8')
end

% first compute within-region averages, for each region
m_f_hp_reg = zeros(R, T);
m_f_pp_reg = zeros(R, T);
m_f_ic_reg = zeros(R, T);
m_f_hs_reg = zeros(R, T);
for r = 1:R
     v_ind = find(v_reg == r); % (Nr,1)
     mat = m_hp(v_ind, :); m_f_hp_reg(r, :) = mean(mat, 1); % (1,T)
     mat = m_pp(v_ind, :); m_f_pp_reg(r, :) = mean(mat, 1);
     mat = m_ic(v_ind, :); m_f_ic_reg(r, :) = mean(mat, 1);
     mat = m_hs(v_ind, :); m_f_hs_reg(r, :) = mean(mat, 1);
end

% intercept
v_ones = ones(T, 1);

% initialise
m_hp_dfac = zeros(N, T);
m_pp_dfac = zeros(N, T);
m_ic_dfac = zeros(N, T);
for ii = 1:N
     % select region corresponding to MSA ii
     r = v_reg(ii);
     v_f_hp_reg = m_f_hp_reg(r, :)'; %(T,1)
     v_f_pp_reg = m_f_pp_reg(r, :)';
     v_f_ic_reg = m_f_ic_reg(r, :)';
     v_f_hs_reg = m_f_hs_reg(r, :)';

     m_f_hp(:, 2) = v_f_hp_reg;
     m_f_pp(:, 2) = v_f_pp_reg;
     m_f_ic(:, 2) = v_f_ic_reg;
     m_f_hs(:, 2) = v_f_hs_reg;

     m_Xi = [v_ones m_quarterly_dummies m_f_hp m_f_pp m_f_ic m_f_hs];
     v_yi = m_hp(ii, :)'; m_hp_dfac(ii, :) = fn_defactoring_cs(v_yi, m_Xi);
     v_yi = m_pp(ii, :)'; m_pp_dfac(ii, :) = fn_defactoring_cs(v_yi, m_Xi);
     v_yi = m_ic(ii, :)'; m_ic_dfac(ii, :) = fn_defactoring_cs(v_yi, m_Xi);

end

% wide2long
v_hp_dfac = m_hp_dfac(:);
v_pp_dfac = m_pp_dfac(:);
v_ic_dfac = m_ic_dfac(:);

% -------------------------------------------------------------------------- %
save(sprintf('%sdata_main_defactored_wide_W%03d.mat', foldername, miles), 'v_id', 'v_time', 'v_reg', 'm_hp_dfac', 'm_pp_dfac', 'm_ic_dfac');

if miles ~= 0 && CD
     % save m_hp (after defactoring) for CD- and alpha-test (currently available only in GAUSS)
     writematrix(m_hp_dfac, sprintf('%shp_defact.txt', foldername), 'Delimiter', 'tab');
end

% -------------------------------------------------------------------------- %
function v_res = fn_defactoring_cs(v_y, m_X)
m_XX = m_X' * m_X;
v_Xy = m_X' * v_y;
v_beta = mldivide(m_XX, v_Xy);
v_res = v_y - (m_X * v_beta);
end
