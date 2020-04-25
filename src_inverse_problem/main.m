%---------------------------%

% Display area D

tn = 0:0.01:2 * pi;

x1 = x_vector(1, tn);
x2 = x_vector(2, tn);
y1 = y_vector(1, tn);
y2 = y_vector(2, tn);

patch([x1, x2],...
      [y1, y2],...
      [.7 .7 .7])
  
global s;

global x2;
global y2;

global r_1;
global r_2;

q_m = 1;

r_1 = cos(s);
r_2 = sin(s);

while abs(norm(q_m)) >= 1e-3
    matrix_initialization

    %---------------------------%

    % Read data from file

    matrix_file = 'A_matrix.txt';
    vector_file = 'y_vector.txt';

    A = dlmread(matrix_file);
    y = dlmread(vector_file);

    %---------------------------%

    % Set exact solution

    m = (size(y, 1) - 3) / 8;

    s11 = -3/2;
    s12 = 3/2;
    s21 = 2;
    s22 = 2;
    s31 = 3;
    s32 = -3;
    s41 = -5/2;
    s42 = -5/2;

    % Exact solution
    % U = @(x,y) (10 * r(x,y,s11,s12)^2 * log(r(x,y,s11,s12)) + 4 * r(x,y,s21,s22)^2 * log(r(x,y,s21,s22)) + ...
    %     3 * r(x,y,s21,s22)^2 * log(r(x,y,s21,s22)) + 2 * r(x,y,s41,s42)^3 * log(r(x,y,s41,s42))) / (8 * pi);

    %---------------------------%

    % Solve system of the equations

    x = A \ y;

    a0 = x(8 * m + 1);
    a1 = x(8 * m + 2);
    a2 = x(8 * m + 3);

    fi_func_1 = zeros(2 * m, 1);
    fi_func_2 = zeros(2 * m, 1);
    psi_func_1 = zeros(2 * m, 1);
    psi_func_2 = zeros(2 * m, 1);

    for ii = 1:2 * m 
        fi_func_1(ii) = x(ii);
        fi_func_2(ii) = x(ii + 2 * m);

        psi_func_1(ii) = x(ii + 4 * m);
        psi_func_2(ii) = x(ii + 6 * m);
    end

    %---------------------------%

    % Find approximate solution U

    U_approx = 0;
    f = 0;

    x = 2;
    y = 0;

    y_1 = 0;
    y_2 = 0;

    r = sqrt((x - y_1)^2 + (y - y_2)^2);

    U = @(x, y) (x^2 + y^2) / 2; %r^2 * log(r) / (8 * pi);

    r_1 = sqrt((x_vector_inner(1, s) - x).^2 + (y_vector_inner(1, s) - y).^2);
    r_2 = sqrt((x_vector_inner(2, s) - x).^2 + (y_vector_inner(2, s) - y).^2);

    H11 =  r_1.^2 .* log(r_1) ./ 4;
    H21 =  -(n1_inner(1, s) .* (x - x_vector_inner(1, s)) .* (1 + 2 .* log(r_1)) +...
             n2_inner(1, s) .* (y - y_vector_inner(1, s)) .* (1 + 2 .* log(r_1))) ./ 4;

    H12 =  r_2.^2 .* log(r_2) ./ 4;
    H22 =  -(n1_inner(2, s) .* (x - x_vector_inner(2, s)) .* (1 + 2 .* log(r_2)) +...
             n2_inner(2, s) .* (y - y_vector_inner(2, s)) .* (1 + 2 .* log(r_2))) ./ 4;

    for ii = 1:2 * m 
        f = f + fi_func_1(ii) * H11(ii) + psi_func_1(ii) * H21(ii) +...
                fi_func_2(ii) * H12(ii) + psi_func_2(ii) * H22(ii);

        U_approx = U_approx + fi_func_1(ii) *...
              H11(ii) + psi_func_1(ii) * H21(ii) +...
                   fi_func_2(ii) * H12(ii) + psi_func_2(ii) * H22(ii);
    end

    f = f / (2 * m) + a0 + a1 * x + a2 * y;

    U_approx = U_approx / (2 * m) + a0 + a1 * x + a2 * y;
    U_exact = U(x, y);

    error = abs(U_approx - U_exact);

    %---------------------------%

    % Display results

    disp([newline char(9) 'RESULTS'])
    disp(['m = ' num2str(m)])
    disp(['f = ' num2str(f)])
    %disp([newline 'U_exact('  num2str(x) ', ' num2str(y) ') = ' num2str(U_exact)])
    %disp(['U_approx(' num2str(x) ', ' num2str(y) ') = ' num2str(U_approx)])
    %disp(['error = ' num2str(error)])

    %---------------------------%

    n = m - 1;

    % Right vector initialization
    b = zeros(2 * m, 1);

    for ii = 1:2 * m

        b_i = 0;

        for k = 1:2 * m
            b_i = b_i + H1(2, 1, ii, k) .* fi_func_1(k) + H2(2, 1, ii, k) .* psi_func_1(k) +...
                (2 * m * R(abs(ii - k), n) .* H1_1(2, ii, k) + H1_2(2, ii, k)) .* fi_func_2(k) +...
                (2 * m * R(abs(ii - k), n) .* H2_1(2, ii, k) + H2_2(2, ii, k)) .* psi_func_2(k);        
        end

        b(ii) = f - (a0 + a1 * x2(ii) + a2 * y2(ii)) - b_i / (2 * m);
    end

    % Matrix initialization
    M = zeros(2 * m, 2 * n + 1);

    for ii = 1:2 * m
        for jj = 1:2 * n + 1

            m_ij = 0;

            for k = 1:2 * m
                m_ij = m_ij + l(jj, s(k), n) .* (N_r_1(ii, jj) .* fi_func_1(k)) +...
                                                 N_r_2(ii, jj) .* psi_func_1(k);
            end

            M(ii, jj) = m_ij / (8 * m);

        end
    end

    lambda = 0.1;

    MT_M = transpose(M) * M;

    q = (MT_M + lambda * eye(size(MT_M)))^(-1) * transpose(M) * b;
    q = reshape(q, 1, size(q, 1));

    q_m = 0;
    for jj = 1:2 * n + 1
        q_m = q_m + q(jj) * l(jj, s(jj), n);
    end
    
    disp(['q_m = ', num2str(q_m)]);
end