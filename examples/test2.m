% https://link.aps.org/doi/10.1103/PhysRevLett.94.230502
% https://math.preview.excalidraw.com/#json=9uW0zpRYOyGenpEQWuhj4,m4lQRlHTyMPie_IzyTYTwQ

IMAG = 0;

d = 3;
n = 4;

density = .2;

rng(0);
if density < 1
  psi = sprandn(d^n, 1, density);
else
  psi = randn(d^n, 1) + IMAG * 1i * randn(d^n, 1);
end
% psi = full(psi);
% psi(d^(n-1)+6:end) = 0;
psi = psi / norm(psi)

function [W] = makeHouseholder(phi)
  zero = zeros(size(phi));
  zero(1,1) = 1;
  dot = zero'*phi;
  if abs(dot) >= (1 - eps) * norm(phi)
    W = NaN;
    return
  end
  if abs(dot) < eps * norm(phi)
    sgn = 1;
  else
    sgn = dot / abs(dot);
  end
  eta = phi - sqrt(phi'*phi) * sgn * zero;
  W = eye(length(phi)) - (2 / (eta'*eta)) * (eta * eta');
end

function seq = makeClubSequence(d, n)
  if n == 1, seq = "-"; return; end
  seq = [];
  seq_ = makeClubSequence(d, n-1);
  for q=0:d-1
    seq = [seq, seq_.insertBefore(1, string(q))];
  end
  seq = [seq, string(repmat('-', 1, n))];
end

function [ctrl, ctrlVal, targ, V] = singleClubHouseholder(term, psi_j, d)
  term = char(term);
  ctrl = -1;
  ctrlVal = 0;
  targ = find(term == '-', 1);
  if targ ~= numel(term) && targ ~= 1
    ctrl = 1:targ - 1;
    ctrlVal = term(ctrl);
  else
    ctrl_ = term > '0';
    if any(ctrl_)
      ctrl = find(ctrl_, 1, 'last');
      ctrlVal = term(ctrl);
    elseif targ > 1
      ctrl = targ - 1;
      ctrlVal = 0;
    end
  end
  ctrlterm = term(term ~= '-');
  phi = zeros(size(psi_j));
  for k = 0:d-1
    t_ = [ctrlterm, char(string(k)), repmat('0', 1, length(term) - length(ctrlterm) - 1)];
    l = zeros(size(psi_j));
    l(base2dec(t_, d) + 1, 1) = 1;
    r = zeros(size(psi_j));
    r(k+1, 1) = 1;
    phi = phi + l' * psi_j * r;
  end
  V = makeHouseholder(phi(1:d));
end


tic
% cir = ReorderState(psi, d, n);
% psi = cir.apply('R', 'N', n, psi, 0, d);
seq = makeClubSequence(d, n);


cir = qclab.QCircuit(n, 0, d);

onedgates = 1;
twodgates = 0;
ndgates = 0;

psis = zeros(numel(psi), numel(seq) + 1);
psis(:, 1) = psi;

printseq = ['   '];
i = 0;
psi_ = psi;
nnzs = [];
affected = [zeros(numel(psi), 1)];
for term = seq
  [c, cv, t, V] = singleClubHouseholder(term, psi_, d);
  if isnan(V)
    continue
  end
  printseq = [printseq, term];
  i = i+1;
  VGate = qclab.qgates.MatrixGate(t-1, V);
  if c == -1
    fVGate = VGate;
    bVGate = VGate.ctranspose();
    onedgates = onedgates + 1;
  elseif isscalar(c)
    fVGate = qclab.qgates.ControlledGate(VGate, c-1, t-1, str2double(cv));
    bVGate = qclab.qgates.ControlledGate(VGate.ctranspose(), c-1, t-1, str2double(cv));
    twodgates = twodgates + 1;
  else
    cv = arrayfun(@(v) str2double(v), cv);
    fVGate = qclab.qgates.MControlledGate(VGate, c-1, t-1, cv);
    bVGate = qclab.qgates.MControlledGate(VGate.ctranspose(), c-1, t-1, cv);
    ndgates = ndgates + 1;
  end
  mat = fVGate.apply('R', 'N', n, eye(length(psi_)), 0, d);
  phi = psi_;
  psi_ = fVGate.apply('R', 'N', n, psi_, 0, d);
  psi_(abs(psi_) < 1e-6) = 0;
  a = reshape(diag(mat) ~= 1, [], 1) + 0;
  a(psi_ == 0 & phi ~= 0) = -1;


  affected = [affected, a];
  nnzs = [nnzs, nnz(psi_)];
  cir.push_back(bVGate);
  psis(:, i+1) = psi_;
end
psis = psis(:, 1:i+1);

figure;
plot(nnzs)

zeros_ = (abs(psis) > 1e-6) + 0;
figure;
h1 = heatmap(zeros_);
h1.XDisplayLabels = printseq;
figure;
H = heatmap(affected);
H.XDisplayLabels = printseq;

PsiPhase = psi_(1,1);

phase = qclab.qgates.Phase(n-1, real(PsiPhase), imag(PsiPhase));
dPhase = qclab.qgates.qudit.SubspaceGate(phase, [1, 0], n-1);

cir.push_back(dPhase);
toc
fprintf('Start nnz: %d, max nnz: %d, gates: %d\n', nnzs(1), max(nnzs), length(nnzs));
% plot(nnzs);

% res = cir.ctranspose().simulate(repmat('0', 1, n)).states;


% if sum(res - psi) < exp(-6)
%   fprintf('The circuit prepares state psi with %d one and %d two qudit gates\n', onedgates, twodgates);
% end
