function tf = is_islanded(pairs, Nnodes)
if nargin < 2 || isempty(Nnodes), Nnodes = max(pairs(:)); end
if isempty(pairs), tf = true; return; end
G = graph(pairs(:,1), pairs(:,2), [], Nnodes);
c = conncomp(G);
tf = (max(c) > 1);
end
