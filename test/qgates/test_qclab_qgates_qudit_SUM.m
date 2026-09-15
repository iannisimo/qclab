classdef test_qclab_qgates_qudit_SUM < matlab.unittest.TestCase
  % Tests for the qudit generalization of CNOT: SUM|c,t> = |c, (c+t) mod d>.
  methods (Test)

    % Explicit d = 3 matrix, built independently index-by-index, for both
    % control < target and control > target orderings.
    function test_matrix_explicit_d3(test)
      d = 3;

      Gct = qclab.qgates.qudit.SUM( 0, 1 );
      expectedCt = zeros(d^2, d^2);
      for c = 0:d-1
        for t = 0:d-1
          inIdx  = c*d + t;
          outIdx = c*d + mod(c+t, d);
          expectedCt(outIdx+1, inIdx+1) = 1;
        end
      end
      test.verifyEqual( Gct.matrix(d), expectedCt );

      Gtc = qclab.qgates.qudit.SUM( 1, 0 );
      expectedTc = zeros(d^2, d^2);
      for t = 0:d-1
        for c = 0:d-1
          inIdx  = t*d + c;
          outIdx = mod(c+t, d)*d + c;
          expectedTc(outIdx+1, inIdx+1) = 1;
        end
      end
      test.verifyEqual( Gtc.matrix(d), expectedTc );
    end

    % SUM must remain unitary (it's a permutation matrix) for any d.
    function test_matrix_is_unitary(test)
      for d = [2, 3, 4, 5]
        G = qclab.qgates.qudit.SUM( 0, 1 );
        mat = G.matrix(d);
        test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 10*eps );
        test.verifyEqual( mat*mat', eye(d^2), 'AbsTol', 10*eps );
      end
    end

    % At d = 2, SUM(control < target) reduces exactly to the standard
    % CNOT matrix.
    function test_backward_compatibility_d2(test)
      G = qclab.qgates.qudit.SUM( 0, 1 );
      expected = [ 1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0 ];
      test.verifyEqual( G.matrix(2), expected );
    end

    % apply() must agree with left-multiplication by matrix() on an
    % adjacent 2-qudit register.
    function test_apply_matches_matrix_adjacent(test)
      d = 3;
      G = qclab.qgates.qudit.SUM( 0, 1 );
      I = qclab.qId(2, false, d);
      applied = G.apply('L', 'N', 2, I, 0, d);
      test.verifyEqual( applied, G.matrix(d) );
    end

    % apply() on non-adjacent qudits pads with identity on the qudits
    % strictly between control and target, for both orderings.
    function test_apply_nonadjacent(test)
      d = 3;

      Gct = qclab.qgates.qudit.SUM( 0, 2 ); % one qudit in between
      I = qclab.qId(3, false, d);
      applied = Gct.apply('L', 'N', 3, I, 0, d);
      expectedCt = zeros(d^3, d^3);
      for c = 0:d-1
        for mid = 0:d-1
          for t = 0:d-1
            inIdx  = c*d*d + mid*d + t;
            outIdx = c*d*d + mid*d + mod(c+t, d);
            expectedCt(outIdx+1, inIdx+1) = 1;
          end
        end
      end
      test.verifyEqual( applied, expectedCt, 'AbsTol', 10*eps );

      Gtc = qclab.qgates.qudit.SUM( 2, 0 ); % control after target
      appliedTc = Gtc.apply('L', 'N', 3, I, 0, d);
      expectedTc = zeros(d^3, d^3);
      for t = 0:d-1
        for mid = 0:d-1
          for c = 0:d-1
            inIdx  = t*d*d + mid*d + c;
            outIdx = mod(c+t, d)*d*d + mid*d + c;
            expectedTc(outIdx+1, inIdx+1) = 1;
          end
        end
      end
      test.verifyEqual( appliedTc, expectedTc, 'AbsTol', 10*eps );
    end

    % qubit / qubits / setQubits bookkeeping.
    function test_qubits(test)
      G = qclab.qgates.qudit.SUM( 2, 0 );
      test.verifyEqual( G.qubit, int64(0) );
      test.verifyEqual( G.qubits, int64([0, 2]) );
      G.setQubits( [1, 3] );
      test.verifyEqual( G.control, int64(1) );
      test.verifyEqual( G.target, int64(3) );
    end

    % equals() compares direction (control vs target order).
    function test_equals(test)
      G1 = qclab.qgates.qudit.SUM( 0, 1 );
      G2 = qclab.qgates.qudit.SUM( 3, 4 );  % same direction, different qudits
      G3 = qclab.qgates.qudit.SUM( 1, 0 );  % swapped direction
      test.verifyTrue( G1.equals(G2) );
      test.verifyFalse( G1.equals(G3) );
    end

    % Constructor rejects control == target.
    function test_control_target_distinct(test)
      threw = false;
      try
        qclab.qgates.qudit.SUM(1, 1);
      catch
        threw = true;
      end
      test.verifyTrue( threw );
    end

    % Building a generalized qudit Bell state: apply the discrete Fourier
    % transform (qclab's qudit generalization of Hadamard) on qudit 2,
    % then SUM(control = 2, target = 0). Since control and target are
    % non-adjacent (qudit 1 sits in between and stays untouched), this
    % also exercises the non-adjacent apply() path end-to-end through a
    % circuit. The result is (1/sqrt(d)) * sum_j |j>_0 |0>_1 |j>_2.
    function test_qudit_bell_state_via_hadamard_and_sum(test)
      d = 3;
      n = 3;
      circuit = qclab.QCircuit( n, 0, d );
      circuit.push_back( qclab.qgates.Hadamard(2) );
      circuit.push_back( qclab.qgates.qudit.SUM(2, 0) );

      psi0 = zeros(d^n, 1);
      psi0(1) = 1; % |0,0,0>
      psi = circuit.matrix * psi0;

      e0 = zeros(d,1); e0(1) = 1;
      expected = zeros(d^n, 1);
      for j = 0:d-1
        ej = zeros(d,1); ej(j+1) = 1;
        expected = expected + kron( ej, kron( e0, ej ) );
      end
      expected = expected / sqrt(d);

      test.verifyEqual( psi, expected, 'AbsTol', 10*eps );
      test.verifyEqual( norm(psi), 1, 'AbsTol', 10*eps );

      % Every basis state with nonzero amplitude has qudit0 == qudit2
      % and qudit1 == 0: this is the entangling signature of SUM.
      for idx = 0:d^n-1
        q0 = floor(idx / d^2);
        q1 = mod(floor(idx / d), d);
        q2 = mod(idx, d);
        if abs(psi(idx+1)) > 10*eps
          test.verifyEqual( q0, q2 );
          test.verifyEqual( q1, 0 );
        end
      end
    end
  end
end
