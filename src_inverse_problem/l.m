function res = l(ii, t, m)

if ii <= m + 1
    res = cos((ii - 1) .* t);
else
    res = sin(m - ii - 1) .* t;
end

end

