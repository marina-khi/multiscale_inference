clear all

% -------------------------------------------------------------------------- %
%cc%
texfolder = './tabfig/';
texfilenam = 'in_tb_hsar_MG_byregion';
miles = 75;
merge = 1; % merge NE+ME and SW+RM

disp(sprintf('  Miles = %d', miles))

% -------------------------------------------------------------------------- %
% load estimates
load(sprintf('./data2/estimates_W%03d.mat', miles), 'results', 'v_id', 'v_reg')
m_theta = results.m_theta;
clear results

% -------------------------------------------------------------------------- %
% remove alpha and sgmsq swap the order of columns
% Note: originally, estimates in model 4 are ordered as [1-Wy 2-1 3-pp 4-ic 5-Wy1 6-y1 7-sgmsq] (N,7)
m_theta_old = m_theta;
N = size(m_theta_old, 1);
m_theta = NaN(N, 6); %cc% [Wy+Wy1 Wy Wy1 y1 pp ic]

% psi_0i + psi_1i
m_theta(:, 1) = m_theta_old(:, 1) + m_theta_old(:, 5);

% Wy Wy1 y1 pp ic
m_theta(:, 2:end) = m_theta_old(:, [1 5 6 3 4]);

% remove MSAs hitting the bounds
v_psi    = m_theta(:, 1);
v_lambda = m_theta(:, 2);
v_psi_1  = m_theta(:, 3);
v_ind = find(abs(v_psi) > 0.994 | abs(v_lambda) > 0.994 | abs(v_psi_1) > 0.994);
if ~isempty(v_ind)
     m_theta(v_ind, :) = [];
     v_id   (v_ind)    = [];
     v_reg  (v_ind)    = [];
end

% retrieve regions' name
Region_Code = [1:8]';
Region_Name = { ...
'New England'   ; ...
'Mideast'       ; ...
'Great Lakes'   ; ...
'Plains'        ; ...
'Southeast'     ; ...
'Southwest'     ; ...
'Rocky Mountain'; ...
'Far West'      };
r_string = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'};
t_reg_unique = table(Region_Code, Region_Name, r_string);

%%% { This part of the code merge regions.
if merge
     % merge NE (r=2) to ME (r=1)
     v_reg(v_reg == 1) = 2;

     % merge RM (r=6) to FW (r=7)
     v_reg(v_reg == 6) = 7;

     % check
     disp(unique(v_reg))

     % fix names as well
     t_reg_unique.Region_Name{1} = 'New England \& Mideast';
     t_reg_unique.Region_Name{2} = 'New England \& Mideast';
     t_reg_unique.Region_Name{6} = 'Southwest \& Rocky Mountain';
     t_reg_unique.Region_Name{7} = 'Southwest \& Rocky Mountain';

     t_reg_unique.r_string{1} = '1 \& 2';
     t_reg_unique.r_string{2} = '1 \& 2';
     t_reg_unique.r_string{6} = '6 \& 7';
     t_reg_unique.r_string{7} = '6 \& 7';
     disp(t_reg_unique)
end
%%% }

% create MG regions
v_reg_unique = unique(v_reg); % (R,1)
v_mg_regions = [v_reg_unique; 0];
n_mg_regions = length(v_mg_regions);

fid = fopen(sprintf('%s%s_W%03d.tex', texfolder, texfilenam, miles), 'w+t');
for r = v_mg_regions'
     if r > 0
          ind = find(t_reg_unique.Region_Code == r);
          reg = t_reg_unique.Region_Name{ind};
          v_ind = find(v_reg == r);
          r_string = t_reg_unique.r_string{ind};
     elseif r == 0
          N = size(m_theta, 1);
          v_ind = 1:N;
          reg = 'US';
          r_string = '';
     end

     % check
     if isempty(v_ind)
          error('No MSAs')
     end

     m_theta_r = m_theta(v_ind, :);
     N_r = size(m_theta_r, 1);
     fprintf(fid, '%6s & %-27s & %3d', r_string, reg, N_r);

     % compute MG estimates, se, and significance stars
     n_col = size(m_theta, 2);
     v_estimate       = zeros(n_col, 1);
     v_estimate_se    = zeros(n_col, 1);
     v_estimate_stars = zeros(n_col, 1);
     for i_col = 1:n_col
          out = fn_MG(m_theta_r(:, i_col));
          v_estimate      (i_col, 1) = out.MG;
          v_estimate_se   (i_col, 1) = out.MG_se;
          v_estimate_stars(i_col, 1) = out.MG_stars;
     end

     % line with estimates and significance-stars
     for i_col = 1:n_col
          estimate       = v_estimate      (i_col, 1);
          estimate_stars = v_estimate_stars(i_col, 1);

          % print estimate
          fprintf(fid, ' & % -5.3f', estimate);

          % print significance stars
          if     estimate_stars == 3 fprintf(fid, '%-14s', '\sym{\ddagger}');
          elseif estimate_stars == 2 fprintf(fid, '%-14s', '\sym{\dagger}');
          elseif estimate_stars == 1 fprintf(fid, '%-14s', '\sym{\ast}');
          elseif estimate_stars == 0 fprintf(fid, '%-14s', '');
          end
     end
     fprintf(fid, ' %s\n', '\\');
     
     % line with se
     fprintf(fid, '%6s & %-27s & %3s', '', '', '');
     for i_col = 1:n_col
          estimate_se    = v_estimate_se   (i_col, 1);
          fprintf(fid, ' & (%5.3f)%13s', estimate_se, '');
     end
     if r == 8 %cc%
          fprintf(fid, ' %s\n', '\\\hline\hline');
     else
          fprintf(fid, ' %s\n', '\\');
     end
end
fclose(fid);
