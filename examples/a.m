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

function vec = dec2base_(val, d, n)
  vec(n) = 0;
  for i = 0:n-1
    i_ = d^i;
    vec(end-i) = mod(floor(val / i_), d);
  end
end

function val = base2dec_(vec, d, n)
  val = 0;
  for i = 0:n-1
    i_ = d^i;
    val = val + (i_ * vec(end-i));
  end
end

function [cir] = ReorderState(psi, n, d)
  X = @qclab.qgates.PauliX;
  SG = @qclab.qgates.qudit.SubspaceGate;
  MCG = @qclab.qgates.MControlledGate;
  cir = qclab.QCircuit(n, 0, d);

  nnz_ = nnz(psi);

  for ptr = 0:nnz_-1
    ptr_dits = dec2base_(ptr, d, n);
    next = find(psi, ptr+1);
    next = next(end) - 1;
    if ptr == next, continue; end
    next_dits = dec2base_(next, d, n);
    for i=numel(ptr_dits):-1:1;
      cur_ptr = ptr_dits(i);
      cur_next = next_dits(i);
      if cur_ptr == cur_next, continue; end
      target = i-1;
      ctrls = [0:target-1, target+1:numel(ptr_dits)-1];
      ctrlStates = next_dits(ctrls+1);
      subspace = [cur_ptr, cur_next];
      subspace = sort(subspace);
      gate = MCG(SG(X, subspace), ctrls, target, ctrlStates);
      psi = gate.apply('R', 'N', n, psi, 0, d);
      cir.push_back(gate);
      next_dits(i) = ptr_dits(i);
    end
  end
end


rng(0);
adj = randomAdjacencyMatrix(d^(n/2), .7);
psi = reshape(adj, [], 1)
psi_ = psi;

nnz_ = nnz(psi);

