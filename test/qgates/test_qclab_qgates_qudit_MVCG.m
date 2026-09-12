classdef test_qclab_qgates_qudit_MVCG < matlab.unittest.TestCase
  methods (Test)

    % MVCG(control < target) applies a different 1-qudit gate on the
    % target for each of the d possible control values: block-diagonal
    % with the i-th gate's matrix on the i-th d x d block.
    function test_matrix_control_lt_target(test)
      d = 3;
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX(), ...
                qclab.qgates.PauliZ() ];
      G = qclab.qgates.qudit.MVCG( 0, 1, gates );

      mat = G.matrix(d);
      expected = blkdiag( qclab.qgates.Identity.matrix(d), ...
                           qclab.qgates.PauliX.matrix(d), ...
                           qclab.qgates.PauliZ.matrix(d) );
      test.verifyEqual( mat, expected );
      test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 10*eps );
    end

    % MVCG(control > target): same semantics with control/target swapped
    % in qubit order; still unitary and equal to its own block form.
    function test_matrix_control_gt_target_is_unitary(test)
      d = 3;
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX(), ...
                qclab.qgates.PauliZ() ];
      G = qclab.qgates.qudit.MVCG( 1, 0, gates );
      mat = G.matrix(d);
      test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 10*eps );
    end

    % apply() must agree with left-multiplication by matrix() on an
    % adjacent 2-qudit register (mirrors the pattern used by the existing
    % qubit SWAP/iSWAP apply tests).
    function test_apply_matches_matrix_adjacent(test)
      d = 3;
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX(), ...
                qclab.qgates.PauliZ() ];
      G = qclab.qgates.qudit.MVCG( 0, 1, gates );

      I = qclab.qId(2, false, d);
      applied = G.apply('L', 'N', 2, I, 0, d);
      test.verifyEqual( applied, G.matrix(d) );
    end

    % apply() on non-adjacent qudits pads with identity on the qudits
    % strictly between control and target.
    function test_apply_nonadjacent(test)
      d = 3;
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX(), ...
                qclab.qgates.PauliZ() ];
      G = qclab.qgates.qudit.MVCG( 0, 2, gates );  % one qudit in between

      I = qclab.qId(3, false, d);
      applied = G.apply('L', 'N', 3, I, 0, d);
      expected = zeros(d^3, d^3);
      for i = 1:d
        expected = expected + kron( qclab.En(i-1, d, false), ...
          kron( qclab.qId(1,false,d), gates(i).matrix(d) ) );
      end
      test.verifyEqual( applied, expected, 'AbsTol', 10*eps );
    end

    % qubit / qubits / setQubits bookkeeping.
    function test_qubits(test)
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX() ];
      G = qclab.qgates.qudit.MVCG( 3, 1, gates );
      test.verifyEqual( G.qubit, int64(1) );
      test.verifyEqual( G.qubits, int64([1, 3]) );
      G.setQubits( [0, 2] );
      test.verifyEqual( G.qubits, int64([0, 2]) );
    end

    % equals() compares direction (control vs target order) and the
    % per-branch gates.
    function test_equals(test)
      gatesA = [ qclab.qgates.Identity(), qclab.qgates.PauliX() ];
      gatesB = [ qclab.qgates.Identity(), qclab.qgates.PauliX() ];
      gatesC = [ qclab.qgates.Identity(), qclab.qgates.PauliZ() ];

      G1 = qclab.qgates.qudit.MVCG( 0, 1, gatesA );
      G2 = qclab.qgates.qudit.MVCG( 0, 1, gatesB );
      G3 = qclab.qgates.qudit.MVCG( 1, 0, gatesB );  % swapped direction
      G4 = qclab.qgates.qudit.MVCG( 0, 1, gatesC );  % different gates

      test.verifyTrue( G1.equals(G2) );
      test.verifyFalse( G1.equals(G3) );
      test.verifyFalse( G1.equals(G4) );
    end

    % Constructor rejects control == target.
    function test_control_target_distinct(test)
      gates = [ qclab.qgates.Identity(), qclab.qgates.PauliX() ];
      threw = false;
      try
        qclab.qgates.qudit.MVCG( 1, 1, gates );
      catch
        threw = true;
      end
      test.verifyTrue( threw );
    end
  end
end
