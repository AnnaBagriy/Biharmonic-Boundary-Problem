function res = N_r_1(ii, jj)

global s;

global n_12;
global n_22;

res = (cos(s(jj)) * n_12(jj) + sin(s(jj)) * n_22(jj)) * (log(r2(2, 1, ii, jj)) + 3);

%disp(r2(2, 2, ii, jj));

end

