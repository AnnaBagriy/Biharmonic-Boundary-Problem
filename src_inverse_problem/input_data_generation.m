
%---------------------------%

% Initialize system of linear equations

matrix_initialization

%---------------------------%

% Read result data from file

matrix_file = 'A_matrix.txt';
vector_file = 'y_vector.txt';

A = dlmread(matrix_file);
y = dlmread(vector_file);

m = (size(y, 1) - 3) / 8;

%---------------------------%

% Solve system of the equations

x = A \ y;

% Retrieve required data

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

disp(['f(x) (x on Г2) = ', num2str(f)]);

%---------------------------%

% Write to file

f_file = 'f_vector.txt';

fid = fopen(f_file, 'wt');

fprintf(fid, '%d\n', f);

fclose(fid);

%---------------------------%
