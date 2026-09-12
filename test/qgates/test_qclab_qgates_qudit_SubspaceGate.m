classdef test_qclab_qgates_qudit_SubspaceGate < matlab.unittest.TestCase
  methods (Test)

    % Embedding a 1-qubit gate into a 2-level subspace {0,1} of a d-level
    % qudit must leave the other levels untouched (identity) and act as
    % the original gate exactly on the chosen subspace.
    function test_embed_adjacent_subspace(test)
      d = 4;
      X = qclab.qgates.PauliX();
      G = qclab.qgates.qudit.SubspaceGate( X, [0, 1], 0 );

      mat = G.matrix(d);
      expected = eye(d);
      expected([1,2],[1,2]) = [0 1; 1 0];
      test.verifyEqual( mat, expected );

      % unitary
      test.verifyEqual( mat'*mat, eye(d), 'AbsTol', eps );
    end

    % The subspace need not be contiguous: {0,2} of a 4-level qudit.
    function test_embed_noncontiguous_subspace(test)
      d = 4;
      X = qclab.qgates.PauliX();
      G = qclab.qgates.qudit.SubspaceGate( X, [0, 2], 0 );

      mat = G.matrix(d);
      expected = eye(d);
      expected([1,3],[1,3]) = [0 1; 1 0];
      test.verifyEqual( mat, expected );
    end

    % A subspace gate on all d levels reduces to the child gate acting
    % on the full qudit (e.g. embedding the qudit-native Hadamard/Fourier
    % gate on {0,...,d-1} is just that gate).
    function test_full_subspace_is_identity_wrapper(test)
      d = 3;
      H = qclab.qgates.Hadamard();
      G = qclab.qgates.qudit.SubspaceGate( H, [0, 1, 2], 0 );
      test.verifyEqual( G.matrix(d), qclab.qgates.Hadamard.matrix(d), 'AbsTol', 10*eps );
    end

    % nbQubits / qubit bookkeeping of the wrapper mirror the wrapped gate.
    function test_nbQubits_and_qubit(test)
      X = qclab.qgates.PauliX();
      G = qclab.qgates.qudit.SubspaceGate( X, [0, 1], 3 );
      test.verifyEqual( G.nbQubits, int64(1) );
      test.verifyEqual( G.qubit, int64(3) );
    end

    % equals() compares subspace and wrapped gate.
    function test_equals(test)
      X = qclab.qgates.PauliX();
      Z = qclab.qgates.PauliZ();
      G1 = qclab.qgates.qudit.SubspaceGate( X, [0, 1], 0 );
      G2 = qclab.qgates.qudit.SubspaceGate( X, [0, 1], 0 );
      G3 = qclab.qgates.qudit.SubspaceGate( X, [0, 2], 0 );
      G4 = qclab.qgates.qudit.SubspaceGate( Z, [0, 1], 0 );

      test.verifyTrue( G1.equals(G2) );
      test.verifyFalse( G1.equals(G3) );
      test.verifyFalse( G1.equals(G4) );
    end

    % toQASM / toTex are explicitly unsupported for qudit gates.
    function test_export_unsupported(test)
      X = qclab.qgates.PauliX();
      G = qclab.qgates.qudit.SubspaceGate( X, [0, 1], 0 );
      test.verifyTrue( test.throws_( @() G.toQASM(1) ) );
      test.verifyTrue( test.throws_( @() G.toTex(1) ) );
    end
  end

  methods (Access = private)
    function bool = throws_(~, fcn)
      bool = false;
      try
        fcn();
      catch
        bool = true;
      end
    end
  end
end
