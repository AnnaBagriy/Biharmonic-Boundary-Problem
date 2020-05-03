
f_file = 'f_vector.txt';

f = dlmread(f_file);
    
global s;
global q_m;
global q_m_prev;

global x2;
global y2;

disp([newline char(9) 'RESULTS'])
    
% Set any default value to start "while" loop
q_m = 1;
q_m_prev = 0;

epsilon = 1e-3;

while abs(norm(q_m)) >= epsilon
    q_m_prev = q_m;
    
    %---------------------------%
    
    % Solve direct problem
    
    matrix_initialization

    %---------------------------%

    size(tn)
    plot(x_vector(1, tn, size(tn, 2)), y_vector(1, tn, size(tn, 2)));
    hold on;
    
    % Read data from file

    matrix_file = 'A_matrix.txt';
    vector_file = 'y_vector.txt';

    A = dlmread(matrix_file);
    y = dlmread(vector_file);

    m = (size(y, 1) - 3) / 8;

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

    % Find f on Г_2 boundary
    
    f = zeros(1, 2 * m);
    
    for ii = 1:2 * m 
        
        f_i = a0 + a1 * x2(ii) + a2 * y2(ii);
        
        for jj = 1:2 * m
            
            f_i = f_i + fi_func_1(jj) * H1(2, 1, ii, jj) / (2 * m) +...
                        fi_func_2(jj) * (H1_1(2, ii, jj) * R(abs(jj - ii), m) + H1_2(2, ii, jj) / (2 * m)) +...
                        psi_func_1(jj) * H2(2, 1, ii, jj) / (2 * m) +...
                        psi_func_2(jj) * (H2_1(2, ii, jj) * R(abs(jj - ii), m) + H2_2(2, ii, jj) / (2 * m));

        end
        
        f(ii) = f_i;
    end

    %---------------------------%

    % Solve incorrect equation
    
    n = 2;

    %-------------%
    
    % Right vector initialization
    
    b = zeros(2 * m, 1);

    for ii = 1:2 * m

        b_i = 0;

        for k = 1:2 * m
            b_i = b_i + H1(2, 1, ii, k) .* fi_func_1(k) + H2(2, 1, ii, k) .* psi_func_1(k) +...
                (2 * m * R(abs(ii - k), n) .* H1_1(2, ii, k) + H1_2(2, ii, k)) .* fi_func_2(k) +...
                (2 * m * R(abs(ii - k), n) .* H2_1(2, ii, k) + H2_2(2, ii, k)) .* psi_func_2(k);        
        end

        b(ii) = f(ii) - (a0 + a1 * x2(ii) + a2 * y2(ii)) - b_i / (2 * m);
    end

    %-------------%
    
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

    %-------------%
    
    % Regularized least-squares method
    
    lambda = 1e-10;
    
    MT_M = transpose(M) * M;
    q = (MT_M + lambda * eye(size(MT_M)))^(-1) * transpose(M) * b;

    %-------------%
    
    % Find q_m vector to improve radial function r
    
    q_m = zeros(1, 2 * m);

    for ii = 1:2 * m
        
        q_m_k = 0;
        
        for jj = 1:2 * n + 1
            q_m_k = q_m_k + q(jj) * l(jj, s(ii), n);
        end
        
        q_m(ii) = q_m_k;
    end
    
    %-------------%
    
    % Display results
    
    disp([newline 'f_2 = ' num2str(f)])
    disp(['||q_m|| = ' num2str(norm(q_m))]);
    
    %-------------%
    
end