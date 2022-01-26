clear all

%-- calculate the inverse distance matrix and save it as csv-file --%

method = 'inverse distance truncated before';

miles = 75; %cc%

% load distance matrix
m_dist = readmatrix('./data2/spatial_matrix_distance.csv'); % (377,377)
N = size(m_dist, 1);

% inverse distance matrix
%% initialise
m_C = NaN(N, N);

%% 0s on the main diagonal
m_C(1:(N + 1):end) = 0;

%% fill off-diagonal elements
for ii = 1:(N - 1)
     for jj = (ii + 1):N
          d_ij = m_dist(ii, jj);

          switch method
          case {'distance'}
               if d_ij <= miles
                    c_ij = 1;
               else
                    c_ij = 0;
               end
          case {'inverse distance'}
               c_ij = 1 / (d_ij^6); %cc%
          case {'inverse distance truncated before'}
               if d_ij <= miles
                    c_ij = 1 / d_ij;
               else
                    c_ij = 0;
               end
          case {'inverse distance truncated after'}
               error('to be completed')
               c_ij = 1 / d_ij;
               if c_ij < 1000 %cc%
                    c_ij = 0;
               end
          otherwise
               error('')
          end
          m_C(ii, jj) = c_ij;
          m_C(jj, ii) = c_ij;
     end
end

%% row normalise 
v_rowsum_C = sum(m_C, 2);
%%% worry about 0s
if any(v_rowsum_C == 0)
     v_ind = find(v_rowsum_C == 0);
     v_rowsum_C(v_ind) = 1;
end
m_rowsum_C = repmat(v_rowsum_C, [1 N]);
m_W = m_C ./ m_rowsum_C;

% save as cvs (NxN matrix)
writematrix(m_W, './data2/W000.csv');
