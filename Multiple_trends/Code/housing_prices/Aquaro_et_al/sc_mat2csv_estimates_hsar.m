clear all
foldername = './data2/'; %cc%

miles = 75; %cc%

% save results
load(sprintf('%sestimates_W%03d.mat', foldername, miles), ...
     'v_id', ...
     'v_reg', ...
     'results', ...
     'v_time' ... included for completeness, not needed
     )

T = length(v_time);
N = length(v_id);
disp(sprintf('  (N, T) = (%d, %d)', N, T))

% check variance
if sum(results.m_sandwich(:) < 0) > 0
     error('negative variance')
end

%% individual estimates
m_theta = results.m_theta(:, 1:(end - 1)); 
m_theta_var = results.m_sandwich(:, 1:(end - 1));
m_theta_se = realsqrt(m_theta_var);
v_sgmsq = results.m_theta(:, end); 

N = size(m_theta, 1);

% header csv-file
s_Wy0     = 'Wy0,Wy0_se,Wy0_sl';
s_1       = 'intercept,intercept_se,intercept_sl';
s_pp      = 'pp,pp_se,pp_sl';
s_ic      = 'ic,ic_se,ic_sl';
s_Wy1     = 'Wy1,Wy1_se,Wy1_sl';
s_y1      = 'y1,y1_se,y1_sl';
s_Wy0_Wy1 = 'Wy0_Wy1,Wy0_Wy1_se,Wy0_Wy1_sl';

fid = fopen(sprintf('%sestimates_W%03d.csv', foldername, miles), 'w+t');
v_cov_psi0i_psi1i = results.v_cov_psi0i_psi1i; 
fprintf(fid, 'fips,region,%s,%s,%s,%s,%s,%s,%s,sgmsq\n', ...
s_Wy0, ...
s_1  , ...
s_pp , ...
s_ic , ...
s_Wy1, ...
s_y1 , ...
s_Wy0_Wy1);
n_col = size(m_theta, 2);
for ii = 1:N
     % print FIPS code and region id
     fprintf(fid, '%d,%d,', v_id(ii), v_reg(ii));

     for i_col = 1:n_col
          estimate = m_theta(ii, i_col);
          estimate_se = m_theta_se(ii, i_col);
          estimate_atval = abs(estimate / estimate_se);
          estimate_stars = (estimate_atval > 1.645) + ...
                           (estimate_atval > 1.96) + ...
                           (estimate_atval > 2.33);

          % print estimate
          fprintf(fid, '%f,', estimate);

          % print robust standard errors (sandwich formula)
          fprintf(fid, '%f,', estimate_se);

          % print significance
          fprintf(fid, '%d,', estimate_stars);

     end
     % print Wy0+Wy1
     psi0i = m_theta(ii, 1); %cc%
     psi1i = m_theta(ii, 5); %cc%
     estimate = psi0i + psi1i;

     psi0i_var = m_theta_var(ii, 1); %cc%
     psi1i_var = m_theta_var(ii, 5); %cc%
     cov01i = v_cov_psi0i_psi1i(ii);
     if psi0i_var + psi1i_var + (2 * cov01i) < 0
          disp(psi0i_var + psi1i_var + (2 * cov01i))
          disp(ii)
     end
     estimate_se = realsqrt(psi0i_var + psi1i_var + (2 * cov01i));

     estimate_atval = abs(estimate / estimate_se);
     estimate_stars = (estimate_atval > 1.645) + ...
                      (estimate_atval > 1.96) + ...
                      (estimate_atval > 2.33);

     % print estimate
     fprintf(fid, '%f,', estimate);

     % print robust standard errors (sandwich formula)
     fprintf(fid, '%f,', estimate_se);

     % print significance
     fprintf(fid, '%d,', estimate_stars);

     % print sgmsq
     fprintf(fid, '%f\n', v_sgmsq(ii));

end
fclose(fid);
