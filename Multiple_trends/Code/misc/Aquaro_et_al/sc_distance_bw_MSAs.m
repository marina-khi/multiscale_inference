clear all

%-- It calculates the distance between any two MSAs (based on Cynthia's "msa_formW.m") --%

% load spatial coordinates (lat-lon)
df_data = readtable('./data1/yang/msa_coord_4matlab.csv', ...
     'Delimiter', ',', ...
     'ReadVariableNames', true); % (377,3)
N = height(df_data); % 377

% _Note_: MSAs in "msa_coord_4matlab.csv" are sorted differently than in "data_main.csv"!!
v_fips_msa_coord_4matlab = df_data.msacode;

%% load data_main
tab = readtable('./data1/yang/data_main.csv', ...
     'Delimiter', ',', ...
     'ReadVariableNames', true); % long format (NT,7)
v_fips_data_main = tab.msacode(1:N); % (N,1)
if ~all(v_fips_msa_coord_4matlab == v_fips_data_main)
     warning(sprintf('  The way in which MSAs are sorted is different across files!\nFixing it now ...\n'))
end

%% re-arrange order
if 1 %cc%
     v_ind = fn_method(v_fips_msa_coord_4matlab, v_fips_data_main);
     if ~all(v_fips_msa_coord_4matlab(v_ind) == v_fips_data_main)
          error('')
     end
     df_data = df_data(v_ind, :);

end
v_id  = df_data.msacode;
v_lat = df_data.msalat;
v_lng = df_data.msalong;

% Convert all decimal degrees to radians
v_lat = v_lat .* pi ./ 180;
v_lng = v_lng .* pi ./ 180;

% compute distance (in miles)
% Ref: BHP Appendix III Haversine formula
R = 3959; % Earth radius in miles
m_dist = zeros(N, N); % (377,377)
for ii = 1:(N - 1)
    for jj = (ii + 1):N
        latdif = v_lat(jj) - v_lat(ii);
        lngdif = v_lng(jj) - v_lng(ii);
        a_ij_1 = sin(latdif / 2)^2;
        a_ij_2 = cos(v_lat(ii)) * cos(v_lat(jj)) * sin(lngdif / 2)^2;
        a_ij = a_ij_1 + a_ij_2;
        c_ij = 2 * atan2(sqrt(a_ij), sqrt(1 - a_ij));
        d_ij = R * c_ij;
        m_dist(ii, jj) = d_ij;
        m_dist(jj, ii) = d_ij;
    end
end

% save as cvs (NxN matrix)
writematrix(m_dist, './data2/spatial_matrix_distance.csv');

%% save as cvs (tidy format)
%fid = fopen('./data2/msa_distance.csv', 'w+t');
%fprintf(fid, 'ii,jj,distance\n');
%for ii = 1:N
%     for jj = 1:N
%          fprintf(fid, '%d,%d,%f\n', v_id(ii), v_id(jj), m_dist(ii, jj));
%     end
%end
%fclose(fid);

% ---------------------------------------------------------------------------- %
% find v_ind such that v_v1(v_ind) == v_v2
% https://it.mathworks.com/matlabcentral/answers/22926-finding-the-indices-of-the-elements-of-one-array-in-another
function v_ind = fn_method(v1, v2)
     v_ind = arrayfun(@(x)find(v1 == x, 1), v2);
end
