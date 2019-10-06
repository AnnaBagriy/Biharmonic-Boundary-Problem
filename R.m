function res = R(k, m)
% Calculates weight function
% for Nystrom method


jj = 1:m - 1;
%r = cos(jj .* (s(ii) - s(k))) ./ jj;
r = cos(jj .* pi .* k / m) ./ jj;

%res = - (1 + 2 * sum(r) + (cos(m * (s(ii) - s(k)))) / m) / (2 * m);
res = - (0.5 + (-1)^k / (2 * m) + sum(r)) / m;
%res = - (sum(r) + (cos(m * (s(ii) - s(k)))) / (2 * m)) / (m);

end

