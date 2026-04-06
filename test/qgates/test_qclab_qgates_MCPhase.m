classdef test_qclab_qgates_MCPhase < matlab.unittest.TestCase
  methods (Test)
    function test_MCPhase(test)
      mcphase = qclab.qgates.MCPhase( [0,1],2,[1,1],pi/2 );
      test.verifyEqual( mcphase.nbQubits, int64(3) );    % nbQubits
      test.verifyFalse( mcphase.fixed );                 % fixed
      test.verifyTrue( mcphase.controlled );             % controlled
      test.verifyEqual( mcphase.controls, int64([0,1]) ); % control
      test.verifyEqual( mcphase.targets, int64(2) );      % target
      test.verifyEqual( mcphase.controlStates, int64([1,1])); % controlState
      
      % matrix
      phmatrix = qclab.qgates.Phase(0,pi/2).matrix;
      P11 = diag([0 0 0 1]);
      I2  = eye(2);
      I4  = eye(4);
      expected = kron(I2, I4 - P11) + kron(phmatrix, P11);
      test.verifyEqual(mcphase.matrix, expected);
      
      % qubit
      test.verifyEqual( mcphase.qubit, int64(0) );
      
      % qubits
      qubits = mcphase.qubits;
      test.verifyEqual( length(qubits), 3 );
      test.verifyEqual( qubits(1), int64(0) );
      test.verifyEqual( qubits(2), int64(1) );
      test.verifyEqual( qubits(3), int64(2) );
      qnew = [5, 3, 1] ;
      mcphase.setQubits( qnew );
      test.verifyEqual( table(mcphase.qubits()).Var1(1), int64(3) );
      test.verifyEqual( table(mcphase.qubits()).Var1(2), int64(5) );
      test.verifyEqual( table(mcphase.qubits()).Var1(3), int64(1) );
      qnew = [0, 1, 2];
      mcphase.setQubits( qnew );
    
      % toQASM
       [T,out] = evalc('mcphase.toQASM(1)'); % capture output to std::out in T
       test.verifyEqual( out, -1 );
    
      % draw gate
       [out] = mcphase.draw(1, 'N');
       test.verifyEqual( out, 0 );
     
       mcphase.setControlStates( [0,0]  );
       [out] = mcphase.draw(1, 'S');
       test.verifyEqual( out, 0 );
     
       mcphase.setControls([3,1]);
       mcphase.update( pi/3 );
       [out] = mcphase.draw(1, 'L');
       test.verifyEqual( out, 0 );
    
       mcphase.setControlStates([1,0]);
       [out] = mcphase.draw(0, 'N');
       test.verifyTrue( isa(out, 'cell') );
       test.verifySize( out, [9, 1] );
     
       mcphase.setControls([0,1]);
       mcphase.update( 0 );
     
       % TeX
       [out] = mcphase.toTex(1, 'N');
       test.verifyEqual( out, 0 );
     
       mcphase.setControlStates( [0,0]  );
       [out] = mcphase.toTex(1, 'S');
       test.verifyEqual( out, 0 );
     
       mcphase.setControls([3,1]);
       mcphase.update( pi/3 );
       [out] = mcphase.toTex(1, 'L');
       test.verifyEqual( out, 0 );
     
       mcphase.setControlStates([1,0]);
       [out] = mcphase.toTex(0, 'N');
       test.verifyTrue( isa(out, 'cell') );
       test.verifySize( out, [3, 1] );
     
       mcphase.setControls([0,1]);
       mcphase.update( 0 );
     
    % gate
       Phase = qclab.qgates.Phase() ;
       test.verifyTrue( mcphase.gate == Phase );
       test.verifyEqual( mcphase.gate.matrix, Phase.matrix );

     % operators == and ~=
      test.verifyTrue( mcphase ~= Phase );
      test.verifyFalse( mcphase == Phase );
      Phase2 = qclab.qgates.MCPhase([0,1],2,[1,0]) ;
      test.verifyTrue( mcphase == Phase2 );
      test.verifyFalse( mcphase ~= Phase2 );

      % set control, target controlState
      mcphase.setControls([3,4]);
      mcphase.setTargets(5);
      test.verifyTrue( mcphase == Phase2 );
      test.verifyFalse( mcphase ~= Phase2 );
      mcphase.setControls([4,6]);
      mcphase.setTargets(1);
      test.verifyTrue( mcphase ~= Phase2 );
      test.verifyFalse( mcphase == Phase2 );      

      % makeFixed, makeVariable
      mcphase.makeFixed() ;
      test.verifyTrue( mcphase.fixed );
      mcphase.makeVariable() ;
      test.verifyFalse( mcphase.fixed );

      % angle, theta, sin, cos
      angle = qclab.QAngle ;
      test.verifyTrue( mcphase.angle == angle );
      test.verifyEqual( mcphase.theta(), 0 );
      test.verifyEqual( mcphase.cos(), 1 );
      test.verifyEqual( mcphase.sin(), 0 );

      % update(angle)
      angle.update( pi/4 );
      mcphase.update( angle ) ;
      test.verifyEqual( mcphase.theta, pi/4, 'AbsTol', eps );
      test.verifyEqual( mcphase.cos, cos(pi/4), 'AbsTol', eps);
      test.verifyEqual( mcphase.sin, sin(pi/4), 'AbsTol', eps);
      
      % update(theta)
       mcphase.update( pi );
       test.verifyEqual( mcphase.theta, pi, 'AbsTol', eps );
       test.verifyEqual( mcphase.cos, cos(pi), 'AbsTol', eps);
       test.verifyEqual( mcphase.sin, sin(pi), 'AbsTol', eps);
       
      % update(cos, sin)
       mcphase.update( cos(pi/4), sin(pi/4) );
       test.verifyEqual( mcphase.theta, pi/4, 'AbsTol', eps );
       test.verifyEqual( mcphase.cos, cos(pi/4), 'AbsTol', eps);
       test.verifyEqual( mcphase.sin, sin(pi/4), 'AbsTol', eps);
       
      % matrix (non-Id)
       mcphase = qclab.qgates.MCPhase(0, 1, 1, cos(pi/4), sin(pi/4) ) ;
       Phase = qclab.qgates.Phase( 0, cos(pi/4), sin(pi/4) ) ;
       E0 = [1 0; 0 0];
       E1 = [0 0; 0 1];
       mat2 = kron(E0,eye(2)) + kron(E1,Phase.matrix);
       test.verifyEqual( mcphase.matrix, mat2, 'AbsTol', eps );
       mcphase.setControls(2);
       mat2 = kron(eye(2),E0) + kron(Phase.matrix,E1);
       test.verifyEqual( mcphase.matrix, mat2, 'AbsTol', eps );
       mcphase.setControlStates(0);
       mat2 = kron(eye(2),E1) + kron(Phase.matrix,E0);
       test.verifyEqual( mcphase.matrix, mat2, 'AbsTol', eps );
       mcphase.setControls(0);
       mat2 = kron(E1,eye(2)) + kron(E0,Phase.matrix);
       test.verifyEqual( mcphase.matrix, mat2, 'AbsTol', eps );
       
      % ctranspose
       mcphase = qclab.qgates.MCPhase(1, 2, 1, pi/3 ) ;
       mcphase = mcphase';
       Phase = qclab.qgates.Phase( 2, pi/3 ) ;
       E0 = [1 0; 0 0];
       E1 = [0 0; 0 1];
       mat2 = kron(E0,eye(2)) + kron(E1,Phase.matrix');
       test.verifyEqual( mcphase.nbQubits, int64(2) );
       test.verifyEqual( mcphase.controls, int64(1) );
       test.verifyEqual( mcphase.targets, int64(2) );
       test.verifyEqual(mcphase.matrix, mat2, 'AbsTol', eps );

     end



    function test_MCPhase_copy(test)
      mcphase = qclab.qgates.MCPhase(0, 1, 1, cos(pi/4), sin(pi/4) ) ;
      cmcphase = copy(mcphase);

      test.verifyEqual(mcphase.qubits, cmcphase.qubits);
      test.verifyEqual(mcphase.theta, cmcphase.theta);

      cmcphase.update( 1 );
      test.verifyEqual(mcphase.qubits, cmcphase.qubits);
      test.verifyNotEqual(mcphase.theta, cmcphase.theta);

      mcphase.setControls(10);
      test.verifyNotEqual(mcphase.qubits, cmcphase.qubits);
      test.verifyNotEqual(mcphase.theta, cmcphase.theta);
    end

  end
end