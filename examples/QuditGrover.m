% Grover search over N = d^n items with n qudits of even dimension d.
% The marked states get their phase flipped by an MCZ onto an ancilla
% prepared in |d/2>, since Z_d |d/2> = omega_d^(d/2) |d/2> = -|d/2>.

d = 6;            % qudit dimension, must be even
n = 2;            % search qudits
marked = [2, 14]; % marked items in {0, ..., d^n - 1}

k = floor(pi/4 * sqrt(d^n / numel(marked)));

H = @qclab.qgates.Hadamard;
MCZ = @qclab.qgates.MCZ;
% base-d digits of x, qudit 0 being the most significant
digits = @(x) mod(floor(x ./ d.^(n-1:-1:0)), d);

oracle = qclab.QCircuit(n+1, 0, d);
for x = marked
  oracle.push_back(MCZ(0:n-1, n, digits(x)));
end

diffuser = qclab.QCircuit(n+1, 0, d);
for q = 0:n-1, diffuser.push_back(H(q)); end
diffuser.push_back(MCZ(0:n-1, n, zeros(1, n)));
for q = 0:n-1, diffuser.push_back(H(q)'); end

grover = qclab.QCircuit(n+1, 0, d);
for q = 0:n-1, grover.push_back(H(q)); end
grover.push_back(qclab.qgates.qudit.INC(n, d/2));
for i = 1:k
  grover.push_back(oracle);
  grover.push_back(diffuser);
end

psi = grover.simulate(repmat('0', 1, n+1)).states;

% the ancilla stays in |d/2>: read the search register off that slice
p = abs(psi(d/2+1:d:end)).^2;
fprintf('k = %d iterations, P(marked) = %.4f\n', k, sum(p(marked+1)));
