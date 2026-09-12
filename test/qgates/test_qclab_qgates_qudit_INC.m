classdef test_qclab_qgates_qudit_INC < matlab.unittest.TestCase
  methods (Test)

    % INC with sumval = 1 is exactly the shift operator X_d.
    % (INC.matrix always returns a sparse matrix internally; compare
    % values via full().)
    function test_INC_sumval1_is_shift_operator(test)
      d = 5;
      inc = qclab.qgates.qudit.INC( 0, 1 );
      test.verifyEqual( full(inc.matrix(d)), qclab.qgates.PauliX.matrix(d) );
    end

    % General sumval realizes X_d^sumval.
    function test_INC_general_sumval(test)
      d = 5;
      X = qclab.qgates.PauliX.matrix(d);
      for s = 0:d-1
        inc = qclab.qgates.qudit.INC( 0, s );
        test.verifyEqual( full(inc.matrix(d)), X^s );
      end
    end

    % Default constructor: qubit 0, sumval 1, adjustable (not fixed).
    function test_defaults(test)
      inc = qclab.qgates.qudit.INC();
      test.verifyEqual( inc.qubit, int64(0) );
      test.verifyFalse( inc.fixed );
      test.verifyEqual( full(inc.matrix(5)), qclab.qgates.PauliX.matrix(5) );
    end

    % update() changes sumval when the gate is adjustable, and is
    % rejected once the gate is fixed.
    function test_update_and_fixed(test)
      inc = qclab.qgates.qudit.INC( 0, 1, false );
      inc.update( 3 );
      test.verifyEqual( full(inc.matrix(5)), qclab.qgates.PauliX.matrix(5)^3 );

      incFixed = qclab.qgates.qudit.INC( 0, 1, true );
      test.verifyTrue( incFixed.fixed );
      test.verifyError( @() incFixed.update(2), 'MATLAB:assertion:failed' );
    end

    % equals() compares sumval.
    function test_equals(test)
      inc1 = qclab.qgates.qudit.INC( 0, 2 );
      inc2 = qclab.qgates.qudit.INC( 3, 2 );  % different qubit, same sumval
      inc3 = qclab.qgates.qudit.INC( 0, 1 );
      test.verifyTrue( inc1.equals(inc2) );
      test.verifyFalse( inc1.equals(inc3) );
    end

    % Backward compatibility: at d = 2, INC(sumval=1) is PauliX.
    function test_backward_compatibility_d2(test)
      inc = qclab.qgates.qudit.INC( 0, 1 );
      test.verifyEqual( full(inc.matrix(2)), [0 1; 1 0] );
    end
  end
end
