function res = radial(s, q_m, size)
% Radial function

global q_m_prev;

if q_m == 1
    res = ones(1, size);
    return;
end

%res = sqrt(cos(s).^2 + 0.25 * sin(s).^2) + (q_m);
res = q_m_prev + q_m;

end

