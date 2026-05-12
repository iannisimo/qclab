% https://link.aps.org/doi/10.1103/PhysRevLett.94.230502

IMAG = 1;

d = 3;
n = 3;

U = randn(d^n) + IMAG * 1i * randn(d^n);
[U, ~] = qr(U);

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

function [ctrl, ctrlVal, targ, V] = singleClubHouseholder(term, psi_j, d)
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

function [circuit] = ClubHouseholder(psi, j, d, n)
  INC = @qclab.qgates.qudit.INC;
  circuit = qclab.QCircuit(n, 0, d);
  j_ = dec2base_(j, d, n);
  for i = n-1:-1:0
    incVal = mod(d - j_(end-i), d);
    if incVal ~= 0, circuit.push_back(INC(n-i-1, incVal, true)); end
  end
  phi = circuit.apply('R', 'N', n, psi, 0, d);
  seq = makeClubSequence(d, n);

  for term = seq
    [c, cv, t, V] = singleClubHouseholder(term, phi, d);
    VGate = qclab.qgates.MatrixGate(t-1, V);
    if c >= 0
      VGate = qclab.qgates.ControlledGate(VGate, c-1, t-1, str2double(cv));
    end
    phi = VGate.apply('R', 'N', n, phi, 0, d);
    circuit.push_back(VGate);
  end

  % TODO use angle
  % phase = 1/phi(1,1);
  % phaseGate = qclab.qgates.Phase(n-1, real(phase), imag(phase));
  % phaseDGate = qclab.qgates.qudit.SubspaceGate(phaseGate, [1, 0], n-1);
  % circuit.push_back(phaseDGate);

  for i = 0:n-1
    if j_(i+1) ~= 0, circuit.push_back(INC(i, j_(i+1), true)); end
  end
end

function circuit = addControls(subcircuit, controls, controlStates, d, n)
  circuit = qclab.QCircuit(n, 0, d);
  for gate = subcircuit.objects
    if isa(gate, 'qclab.qgates.QControlledGate2')
      gate_ = gate.gate();
      controls_ = [controls, gate.control+1];
      controlStates_ = [controlStates, gate.controlState];
      target = gate.target+1;
    elseif isa(gate, 'qclab.qgates.QMultiControlledGate')
      gate_ = gate.gate();
      controls_ = [controls, gate.controls+1];
      controlStates_ = [controlStates, gate.controlStates];
      target = gate.targets+1;
    else
      gate_ = gate;
      controls_ = controls;
      controlStates_ = controlStates;
      target = gate.qubit+1;
    end

    if isempty(controls_)
      if ismethod(gate_, 'setQubits')
        gate_.setQubits(target);
      else
        gate_.setQubit(target);
      end
      circuit.push_back(gate_);
    elseif isscalar(controls_)
      circuit.push_back(qclab.qgates.ControlledGate(gate_, controls_, target, controlStates_));
    else
      circuit.push_back(qclab.qgates.MControlledGate(gate_, controls_, target, controlStates_));
    end
  end
end

function circuit = triangle(U, d, n)
  circuit = qclab.QCircuit(n, 0, d);
  if n == 1
    for i = 0:d-2
      HS = eye(d);
      phi = U(i+1:end, i+1);
      HS(i+1:end, i+1:end) = makeHouseholder(phi);
      U = HS * U;
      circuit.push_back(qclab.qgates.MatrixGate(0, HS));
    end
    return
  end
  subcircuit = triangle(U(1:d^(n-1),1:d^(n-1)), d, n-1);
  cir = addControls(subcircuit, [], [], d, n);
  U = cir.apply('R', 'N', n, U, 0, d);
  for g = cir.objects
    circuit.push_back(g);
  end
  b_size = floor(d^(n-1)); % block size
  for k=0:d-2
    c_s = floor(b_size * k); % first column of block
    c_e = floor(b_size * (k+1) - 1); % last column of block
    for j = c_s:c_e
      for l = (k+1):d-1 % Block-row
        r_s = floor(b_size * l); % first row of block
        r_e = floor(b_size * (l+1) - 1); % last row of block
        j_ = dec2base_(j, d, n);
        j_(1) = mod(k+l, d);
        nz = base2dec_(j_, d, n);
        subcircuit = ClubHouseholder(U(r_s+1:r_e+1, j+1), nz, d, n-1);
        controls = 0;
        controlStates = mod(l, d);
        cir = addControls(subcircuit, controls, controlStates, d, n);
        U = cir.apply('R', 'N', n, U, 0, d);
        for g = cir.objects
          circuit.push_back(g);
        end
        
      end
      cir = qclab.QCircuit(n, 0, d);
      controls = (1:n-1);
      controlStates = dec2base_(j, d, n);
      controlStates = controlStates(2:end);
      phi = U(j+1:d^(n-1):d^n, j+1);
      HR = eye(d);
      HR(k+1:end, k+1:end) = makeHouseholder(phi);
      gate = qclab.qgates.MatrixGate(0, HR);
      if isscalar(controls)
        cir.push_back(qclab.qgates.ControlledGate(gate, controls, 0, controlStates));
      else
        cir.push_back(qclab.qgates.MControlledGate(gate, controls, 0, controlStates));
      end
      U = cir.apply('R', 'N', n, U, 0, d);
      for g = cir.objects
        circuit.push_back(g);
      end
    end
    subcircuit = triangle(U(b_size * (k+1) + 1:b_size * (k+2), b_size * (k+1) + 1:b_size * (k+2)), d, n-1);
    cir = addControls(subcircuit, 0, k+1, d, n);
    U = cir.apply('R', 'N', n, U, 0, d);
    for g = cir.objects
      circuit.push_back(g);
    end

  end
  
end

cir = triangle(U, d, n);
res = cir.apply('R', 'N', n, U, 0, d);
norm(abs(res) - eye(size(res)))

% TODO missing global phases
