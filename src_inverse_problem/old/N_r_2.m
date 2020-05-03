function res = N_r_2(ii, jj)

global s;

global x1;
global y1;
global x2;
global y2;

res = (cos(s(jj)) * (x2(ii) - x1(jj)) + sin(s(jj)) * (y2(ii) - y1(jj))) * (log(r2(2, 1, ii, jj)) + 1);

end