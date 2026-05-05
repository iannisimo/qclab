% https://link.aps.org/doi/10.1103/PhysRevLett.94.230502


IMAG = 0;

d = 2;
n = 2;

psi = randn(d^n, 1) + IMAG * 1i * randn(d^n, 1);
psi = psi / norm(psi);

function [W] = makeHouseholder(phi)
  zero = zeros(size(phi));
  zero(1,1) = 1;
  eta = phi - sqrt(phi'*phi) * (zero'*phi / abs(zero'*phi)) * zero;
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

function [ctrl, ctrlVal, targ, V] = singleClubHouseholder(term, psi_j, d);
  term = char(term);
  ctrl = -1;
  ctrlVal = 0;
  targ = find(term == '-', 1);
  ctrl_ = term > '0';
  if any(ctrl_)
    ctrl = find(ctrl_, 1, 'last');
    ctrlVal = term(ctrl);
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

seq = makeClubSequence(d, n);

cir = qclab.QCircuit(n, 0, d);

psi_ = psi;
for term = seq
  [c, cv, t, V] = singleClubHouseholder(term, psi_, d);
  VGate = qclab.qgates.MatrixGate(t-1, V);
  if c == -1
    fVGate = VGate;
    bVGate = VGate.ctranspose();
  else
    fVGate = qclab.qgates.ControlledGate(VGate, c-1, t-1, str2num(cv));
    bVGate = qclab.qgates.ControlledGate(VGate.ctranspose(), c-1, t-1, str2num(cv));
  end
  psi_ = fVGate.apply('R', 'N', n, psi_, 0, d);
  cir.push_back(bVGate);
end

%> TODO add phase gate at the end/start

psi
psi_
res = cir.ctranspose().simulate(repmat('0', 1, n)).states

sum(res - psi)
