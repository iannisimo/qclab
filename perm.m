d = 3;
n = 4;

rng(0);
function A = randomAdjacencyMatrix(n, p, directed, weighted)
% RANDOMADJACENCYMATRIX  Generate a random adjacency matrix.
%
%   A = randomAdjacencyMatrix(N) generates an N×N unweighted, undirected
%   random adjacency matrix with edge probability 0.5.
%
%   A = randomAdjacencyMatrix(N, P) uses edge probability P (0 < P <= 1).
%
%   A = randomAdjacencyMatrix(N, P, DIRECTED) if DIRECTED is true, the
%   matrix is not forced to be symmetric. Default: false.
%
%   A = randomAdjacencyMatrix(N, P, DIRECTED, WEIGHTED) if WEIGHTED is
%   true, edge weights are uniform random values in (0, 1]. Default: false.
%
% Inputs:
%   n        - Number of nodes (positive integer)
%   p        - Edge probability, default = 0.5
%   directed - Logical, default = false
%   weighted - Logical, default = false
%
% Output:
%   A        - N×N adjacency matrix (sparse)
%
% Example:
%   A = randomAdjacencyMatrix(10, 0.3);          % undirected, unweighted
%   A = randomAdjacencyMatrix(10, 0.3, true);    % directed
%   A = randomAdjacencyMatrix(10, 0.3, true, true); % directed + weighted

    % --- Defaults ---------------------------------------------------------
    if nargin < 2 || isempty(p),        p        = 0.5;  end
    if nargin < 3 || isempty(directed), directed = false; end
    if nargin < 4 || isempty(weighted), weighted = false; end

    % --- Input validation -------------------------------------------------
    assert(isnumeric(n) && isscalar(n) && n > 0 && floor(n) == n, ...
        'n must be a positive integer.');
    assert(isnumeric(p) && isscalar(p) && p > 0 && p <= 1, ...
        'p must be a scalar in (0, 1].');

    % --- Generate matrix --------------------------------------------------
    % TODO: swap rand for your preferred RNG seed if reproducibility needed
    R = rand(n, n);

    if weighted
        A = R .* (R < p);          % weights in (0, p]; rescale if needed
    else
        A = double(R < p);         % binary edges
    end

    % No self-loops
    A(1:n+1:end) = 0;

    % Force symmetry for undirected graphs
    if ~directed
        A = triu(A, 1);            % keep upper triangle ...
        A = A + A';                % ... mirror to lower
    end
end

% adj = randomAdjacencyMatrix(4, .7)
% vec = reshape(adj, [], 1);
% vec

cir = qclab.QCircuit(n, 0, d);

X = @qclab.qgates.PauliX;
% MC = @qclab.qgates.MControlledGate;

% cir.push_back(MC(X, [0,1], 2, [1,0]));
% cir.push_back(X(2));
cir.push_back(MC(X, [n-1], 0, [1]));
cir.push_back(MC(X, [n-1], 0, [2]));
cir.push_back(MC(X, [n-1], 0, [2]));
a = cir.matrix();
cir.push_back(MC(X, [0], 1, [2]));

cir.draw;

heatmap(a - cir.matrix())
