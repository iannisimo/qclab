classdef test_qclab_QCircuit_qudit < matlab.unittest.TestCase
  % Tests for the qudit generalization of qclab.QCircuit: the dimension
  % parameter d, the isQudit()/canInsert() guard, and backward
  % compatibility of the default d = 2 behavior.
  methods (Test)

    % Default dimension is 2, and is unaffected by the addition of the
    % new constructor argument (backward compatibility).
    function test_default_dimension_is_2(test)
      circuit = qclab.QCircuit( 2 );
      circuit.push_back( qclab.qgates.PauliX(0) );
      circuit.push_back( qclab.qgates.Hadamard(1) );
      expected = kron( qclab.qgates.PauliX.matrix(2), ...
                        qclab.qgates.Hadamard.matrix(2) );
      test.verifyEqual( circuit.matrix, expected, 'AbsTol', 10*eps );
    end

    % A circuit built with d = 3 produces a 3^nbQubits-sized unitary, and
    % single-qudit gates on disjoint qudits combine via the Kronecker
    % product exactly as in the qubit case.
    function test_qutrit_circuit_matrix(test)
      d = 3;
      circuit = qclab.QCircuit( 2, 0, d );
      circuit.push_back( qclab.qgates.PauliX(0) );
      circuit.push_back( qclab.qgates.Hadamard(1) );

      mat = circuit.matrix;
      test.verifySize( mat, [d^2, d^2] );
      test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 10*eps );

      expected = kron( qclab.qgates.PauliX.matrix(d), ...
                        qclab.qgates.Hadamard.matrix(d) );
      test.verifyEqual( mat, expected, 'AbsTol', 10*eps );
    end

    % canInsert() / push_back() must reject a gate whose matrix was not
    % generalized to arbitrary d (e.g. RotationX, which only received the
    % cosmetic apply(obj, ~) signature change) as soon as d > 2, while
    % still accepting it for d = 2.
    function test_canInsert_rejects_non_qudit_gate_for_d_gt_2(test)
      rx = qclab.qgates.RotationX( 0, pi/4 );

      circuit2 = qclab.QCircuit( 1, 0, 2 );
      test.verifyTrue( circuit2.canInsert( rx ) );
      circuit2.push_back( rx );  % must not throw

      circuit3 = qclab.QCircuit( 1, 0, 3 );
      test.verifyFalse( circuit3.canInsert( rx ) );
      test.verifyError( @() circuit3.push_back( rx ), 'MATLAB:assertion:failed' );
    end

    % canInsert() / push_back() must accept a gate that was actually
    % generalized (PauliX, Hadamard, ...) into a d > 2 circuit.
    function test_canInsert_accepts_qudit_gate_for_d_gt_2(test)
      x = qclab.qgates.PauliX( 0 );
      circuit = qclab.QCircuit( 1, 0, 4 );
      test.verifyTrue( circuit.canInsert( x ) );
      circuit.push_back( x );  % must not throw
      test.verifyEqual( circuit.matrix, qclab.qgates.PauliX.matrix(4) );
    end

    % The new qudit-native gate classes (SubspaceGate, INC, MVCG) must
    % also be accepted into a d > 2 circuit.
    function test_canInsert_accepts_new_qudit_gate_classes(test)
      d = 4;
      circuit = qclab.QCircuit( 2, 0, d );

      subspace = qclab.qgates.qudit.SubspaceGate( qclab.qgates.PauliX(), [0,1], 0 );
      inc = qclab.qgates.qudit.INC( 1, 1 );
      test.verifyTrue( circuit.canInsert( subspace ) );
      test.verifyTrue( circuit.canInsert( inc ) );

      circuit.push_back( subspace );
      circuit.push_back( inc );
      mat = circuit.matrix;
      test.verifySize( mat, [d^2, d^2] );
      test.verifyEqual( mat'*mat, eye(d^2), 'AbsTol', 10*eps );
    end

    % canInsert() rejects a gate placed on a qudit index beyond the
    % circuit's own qudit count, independently of d.
    function test_canInsert_rejects_out_of_range_qubit(test)
      x = qclab.qgates.PauliX( 5 );
      circuit = qclab.QCircuit( 2, 0, 3 );
      test.verifyFalse( circuit.canInsert( x ) );
    end
  end
end
