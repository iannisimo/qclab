classdef test_qclab_qgates_qudit_generalization < matlab.unittest.TestCase
  % Regression + correctness tests for the d-dimensional generalization of
  % the existing qubit gates (Identity, PauliX, PauliY, PauliZ, Hadamard,
  % Phase) introduced on the qudit branch. See Chapter "Extending QCLAB to
  % Qudits" for the closed-form definitions being checked here.
  methods (Test)

    % Every generalized gate must reduce exactly to its original qubit
    % matrix at d = 2, both when called with an explicit d = 2 and with no
    % argument at all (the two call conventions used throughout the
    % existing qubit-only test suite).
    function test_backward_compatibility_d2(test)
      test.verifyEqual( qclab.qgates.Identity.matrix(2), eye(2) );
      test.verifyEqual( qclab.qgates.Identity.matrix,    eye(2) );

      test.verifyEqual( qclab.qgates.PauliX.matrix(2), [0 1; 1 0] );
      test.verifyEqual( qclab.qgates.PauliX.matrix,    [0 1; 1 0] );

      test.verifyEqual( qclab.qgates.PauliY.matrix(2), [0 -1i; 1i 0] );
      test.verifyEqual( qclab.qgates.PauliY.matrix,    [0 -1i; 1i 0] );

      test.verifyEqual( qclab.qgates.PauliZ.matrix(2), [1 0; 0 -1] );
      test.verifyEqual( qclab.qgates.PauliZ.matrix,    [1 0; 0 -1] );

      sqrt2 = 1/sqrt(2);
      Href = [sqrt2, sqrt2; sqrt2, -sqrt2];
      test.verifyEqual( qclab.qgates.Hadamard.matrix(2), Href, 'AbsTol', eps );
      test.verifyEqual( qclab.qgates.Hadamard.matrix,    Href, 'AbsTol', eps );

      P = qclab.qgates.Phase(0, pi/3);
      test.verifyEqual( P.matrix(2), [1, 0; 0, exp(1i*pi/3)], 'AbsTol', eps );
      test.verifyEqual( P.matrix,    [1, 0; 0, exp(1i*pi/3)], 'AbsTol', eps );
    end

    % Every generalized gate must remain unitary for d > 2.
    function test_unitary_for_d_greater_than_2(test)
      for d = [3, 4, 5, 7]
        test.verifyUnitary_( qclab.qgates.Identity.matrix(d), d );
        test.verifyUnitary_( qclab.qgates.PauliX.matrix(d), d );
        test.verifyUnitary_( qclab.qgates.PauliY.matrix(d), d );
        test.verifyUnitary_( qclab.qgates.PauliZ.matrix(d), d );
        test.verifyUnitary_( qclab.qgates.Hadamard.matrix(d), d );
        P = qclab.qgates.Phase(0, 0.7);
        test.verifyUnitary_( P.matrix(d), d );
      end
    end

    % PauliX generalizes to the cyclic shift operator |k> -> |k+1 mod d>.
    function test_PauliX_is_shift_operator(test)
      d = 5;
      X = qclab.qgates.PauliX.matrix(d);
      test.verifyEqual( X, circshift(eye(d), 1, 1) );
      % acting on basis vector e_k gives e_{k+1 mod d}
      for k = 0:d-1
        ek = zeros(d,1); ek(k+1) = 1;
        ekp1 = zeros(d,1); ekp1(mod(k+1,d)+1) = 1;
        test.verifyEqual( X*ek, ekp1 );
      end
    end

    % PauliZ generalizes to the clock operator diag(exp(2pi i k / d)).
    function test_PauliZ_is_clock_operator(test)
      d = 5;
      Z = qclab.qgates.PauliZ.matrix(d);
      k = (0:d-1)';
      test.verifyEqual( Z, diag(exp(2i*pi*k/d)), 'AbsTol', eps );
    end

    % PauliY is defined as i * X_d * Z_d, consistent with the shift/clock
    % (Weyl-Heisenberg) generalization.
    function test_PauliY_consistent_with_X_and_Z(test)
      for d = [3, 4, 5]
        X = qclab.qgates.PauliX.matrix(d);
        Z = qclab.qgates.PauliZ.matrix(d);
        Y = qclab.qgates.PauliY.matrix(d);
        test.verifyEqual( Y, 1i * X * Z, 'AbsTol', 10*eps );
      end
    end

    % For d > 2 the "Hadamard" class computes the discrete Fourier
    % transform F_d = exp(2*pi*i*n'*n/d)/sqrt(d).
    function test_Hadamard_is_discrete_fourier_transform(test)
      for d = [3, 4, 5]
        H = qclab.qgates.Hadamard.matrix(d);
        n = (0:d-1)';
        Fd = exp(2i*pi*(n*n')/d) / sqrt(d);
        test.verifyEqual( H, Fd, 'AbsTol', 10*eps );
      end
    end

    % Phase generalizes to a phase on the top energy level only:
    % diag(1, ..., 1, exp(i*theta)).
    function test_Phase_applies_to_top_level_only(test)
      theta = 1.23456;
      P = qclab.qgates.Phase(0, theta);
      for d = [3, 4, 5]
        mat = P.matrix(d);
        expected = eye(d);
        expected(d,d) = exp(1i*theta);
        test.verifyEqual( mat, expected, 'AbsTol', eps );
      end
    end
  end

  methods (Access = private)
    function verifyUnitary_(test, mat, d)
      test.verifyEqual( mat'*mat, eye(d), 'AbsTol', 10*eps );
      test.verifyEqual( mat*mat', eye(d), 'AbsTol', 10*eps );
    end
  end
end
