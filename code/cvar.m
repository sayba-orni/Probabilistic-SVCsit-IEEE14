function c = cvar(x, alpha)
x = sort(x(:)); k = max(1, ceil(alpha*numel(x)));
c = mean(x(k:end));
end
