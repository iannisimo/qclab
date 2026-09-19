classdef test_qclab_qgates_qudit_SWAP < matlab.unittest.TestCase
  % Tests for the qudit generalization of SWAP: SWAP|q0,q1> = |q1,q0>.
  methods (Test)

    % Explicit d = 3 matrix, built independently index-by-index.
    function test_matrix_explicit_d3(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 1 );

      expected = zeros(d^2, d^2);
      for q0 = 0:d-1
        for q1 = 0:d-1
          inIdx  = q0*d + q1;
          outIdx = q1*d + q0;
          expected(outIdx+1, inIdx+1) = 1;
        end
      end
      test.verifyEqual( G.matrix(d), expected );
    end

    % At d = 2, matrix() reduces exactly to the standard SWAP matrix.
    function test_backward_compatibility_d2(test)
      G = qclab.qgates.SWAP( 0, 1 );
      expected = [ 1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1 ];
      test.verifyEqual( G.matrix(2), expected );
      test.verifyEqual( G.matrix, expected ); % default d = 2
    end

    % SWAP is a permutation matrix, hence unitary, for any d.
    function test_matrix_is_unitary(test)
      G = qclab.qgates.SWAP( 0, 1 );
      for d = [2, 3, 4, 5]
        mat = G.matrix(d);
        test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 1e-10 );
        test.verifyEqual( mat*mat', eye(d^2), 'AbsTol', 1e-10 );
      end
    end

    % SWAP is symmetric and its own inverse for any d.
    function test_matrix_symmetric_involution(test)
      G = qclab.qgates.SWAP( 0, 1 );
      for d = [2, 3, 4, 5]
        mat = G.matrix(d);
        test.verifyEqual( mat, mat.' );
        test.verifyEqual( mat*mat, eye(d^2) );
      end
    end

    % fixed / controlled flags, and isQudit (checked indirectly through
    % QCircuit.canInsert, since isQudit is a protected method).
    function test_flags(test)
      G = qclab.qgates.SWAP( 0, 1 );
      test.verifyTrue( G.fixed );
      test.verifyFalse( G.controlled );

      circuit = qclab.QCircuit( 2, 0, 3 ); % d = 3 qudit circuit
      test.verifyTrue( circuit.canInsert( G ) );
    end

    % apply() must agree with left-multiplication by matrix() on an
    % adjacent 2-qudit register.
    function test_apply_matches_matrix_adjacent(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 1 );
      I = qclab.qId(2, false, d);
      applied = G.apply('L', 'N', 2, I, 0, d);
      test.verifyEqual( applied, G.matrix(d), 'AbsTol', 1e-10 );
    end

    % apply() from the right also matches matrix() (SWAP is symmetric).
    function test_apply_matches_matrix_right(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 1 );
      I = qclab.qId(2, false, d);
      applied = G.apply('R', 'N', 2, I, 0, d);
      test.verifyEqual( applied, G.matrix(d), 'AbsTol', 1e-10 );
    end

    % Transposed / conjugate-transposed application: since the qudit SWAP
    % matrix is real and symmetric, 'T' and 'C' must give the same result
    % as 'N'.
    function test_apply_op_transpose_conjugate(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 1 );
      I = qclab.qId(2, false, d);
      appliedN = G.apply('L', 'N', 2, I, 0, d);
      appliedT = G.apply('L', 'T', 2, I, 0, d);
      appliedC = G.apply('L', 'C', 2, I, 0, d);
      test.verifyEqual( appliedT, appliedN, 'AbsTol', 1e-10 );
      test.verifyEqual( appliedC, appliedN, 'AbsTol', 1e-10 );
    end

    % apply() on non-adjacent qudits pads with identity on the qudit(s)
    % strictly between the two swapped qudits.
    function test_apply_nonadjacent(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 2 ); % one qudit in between
      I = qclab.qId(3, false, d);
      applied = G.apply('L', 'N', 3, I, 0, d);

      expected = zeros(d^3, d^3);
      for q0 = 0:d-1
        for mid = 0:d-1
          for q2 = 0:d-1
            inIdx  = q0*d*d + mid*d + q2;
            outIdx = q2*d*d + mid*d + q0;
            expected(outIdx+1, inIdx+1) = 1;
          end
        end
      end
      test.verifyEqual( applied, expected, 'AbsTol', 1e-10 );
    end

    % apply() offset parameter shifts which qudits are acted on.
    function test_apply_offset(test)
      d = 3;
      G = qclab.qgates.SWAP( 0, 1 );
      I = qclab.qId(3, false, d);
      applied = G.apply('L', 'N', 3, I, 1, d);
      expected = kron( qclab.qId(1, false, d), G.matrix(d) );
      test.verifyEqual( applied, expected, 'AbsTol', 1e-10 );
    end

    % Round-trip through a qudit circuit: applying SWAP twice on the same
    % pair of qudits is the identity, for any d.
    function test_double_swap_is_identity_via_circuit(test)
      d = 3;
      n = 2;
      circuit = qclab.QCircuit( n, 0, d );
      circuit.push_back( qclab.qgates.SWAP(0, 1) );
      circuit.push_back( qclab.qgates.SWAP(0, 1) );
      test.verifyEqual( circuit.matrix, eye(d^n), 'AbsTol', 1e-10 );
    end

    % SWAP correctly exchanges the state of two distinct qudits.
    function test_swap_exchanges_basis_states(test)
      d = 4;
      n = 2;
      circuit = qclab.QCircuit( n, 0, d );
      circuit.push_back( qclab.qgates.SWAP(0, 1) );

      a = 1; b = 3; % |a>_0 |b>_1 -> |b>_0 |a>_1
      psi0 = zeros(d^n, 1);
      psi0( a*d + b + 1 ) = 1;
      psi = circuit.matrix * psi0;

      expected = zeros(d^n, 1);
      expected( b*d + a + 1 ) = 1;
      test.verifyEqual( psi, expected, 'AbsTol', 1e-10 );
    end

  end
end
