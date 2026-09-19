classdef test_qclab_qgates_qudit_GCG < matlab.unittest.TestCase
  % Tests for the qudit Generator Controlled Gate GCG|c,t> = G^c applied
  % to the target qudit t.
  methods (Test)

    % GCG with generator = PauliX (the qudit cyclic shift) must coincide
    % exactly with SUM, which is defined as its generalization:
    % SUM|c,t> = |c, (c+t) mod d>.
    function test_GCG_pauliX_equals_SUM_matrix(test)
      for d = [2, 3, 4, 5]
        Gct = qclab.qgates.qudit.GCG( 0, 1, qclab.qgates.PauliX() );
        Sct = qclab.qgates.qudit.SUM( 0, 1 );
        test.verifyEqual( Gct.matrix(d), Sct.matrix(d), 'AbsTol', 10*eps );

        Gtc = qclab.qgates.qudit.GCG( 1, 0, qclab.qgates.PauliX() );
        Stc = qclab.qgates.qudit.SUM( 1, 0 );
        test.verifyEqual( Gtc.matrix(d), Stc.matrix(d), 'AbsTol', 10*eps );
      end
    end

    % Same equivalence, but through apply() on an adjacent 2-qudit
    % register, for both control/target orderings.
    function test_GCG_pauliX_equals_SUM_apply(test)
      d = 3;
      I = qclab.qId(2, false, d);

      Gct = qclab.qgates.qudit.GCG( 0, 1, qclab.qgates.PauliX() );
      Sct = qclab.qgates.qudit.SUM( 0, 1 );
      test.verifyEqual( Gct.apply('L', 'N', 2, I, 0, d), ...
                         Sct.apply('L', 'N', 2, I, 0, d), 'AbsTol', 10*eps );

      Gtc = qclab.qgates.qudit.GCG( 1, 0, qclab.qgates.PauliX() );
      Stc = qclab.qgates.qudit.SUM( 1, 0 );
      test.verifyEqual( Gtc.apply('L', 'N', 2, I, 0, d), ...
                         Stc.apply('L', 'N', 2, I, 0, d), 'AbsTol', 10*eps );
    end

    % Same equivalence, non-adjacent qudits (exercises GCG's general
    % non-nearest-neighbor apply path).
    function test_GCG_pauliX_equals_SUM_apply_nonadjacent(test)
      d = 3;
      I = qclab.qId(3, false, d);

      G = qclab.qgates.qudit.GCG( 0, 2, qclab.qgates.PauliX() );
      S = qclab.qgates.qudit.SUM( 0, 2 );
      test.verifyEqual( G.apply('L', 'N', 3, I, 0, d), ...
                         S.apply('L', 'N', 3, I, 0, d), 'AbsTol', 10*eps );
    end

    % At d = 2, GCG(PauliX) also reduces to the standard CNOT matrix,
    % consistent with SUM's own d = 2 backward compatibility.
    function test_GCG_pauliX_backward_compatibility_d2(test)
      G = qclab.qgates.qudit.GCG( 0, 1, qclab.qgates.PauliX() );
      expected = [ 1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0 ];
      test.verifyEqual( G.matrix(2), expected );
    end

  end
end
