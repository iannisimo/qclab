classdef test_qclab_HandleCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_HandleCircuit(test)
      X = qclab.qgates.PauliX ;
      Circ1 = qclab.QCircuit(1) ;
      Circ1.push_back(X) ;
      H = qclab.HandleCircuit(Circ1) ;

      test.verifyEqual( H.nbQubits, int64(1) );     % nbQubits
      test.verifyFalse( H.fixed );                   % fixed
      test.verifyFalse( H.controlled );             % controlled
      test.verifyEqual( H.qubit, int64(0) );        % qubit
      test.verifyEqual( H.offset, int64(0) );       % offset
      test.verifySameHandle( H.circuitHandle, Circ1 );     % circuitHandle

      % qubits
      qubits = H.qubits;
      test.verifyEqual( length(qubits), 1 );
      test.verifyEqual( qubits(1), int64(0) );

      % matrix
      test.verifyEqual( H.matrix, X.matrix );

      % toQASM
      [T,out] = evalc('H.toQASM(1)'); % capture output to std::out in T
      test.verifyEqual( out, 0 );
      QASMstring = sprintf(['\n\nQASM output:\n\n' ...
                      'OPENQASM 2.0;\n' ...
                      'include "qelib1.inc";\n\n' ...
                      'qreg q[1];\n' ...
                      'x q[0];\n' ...
                    ]);
      test.verifyEqual(T(1:length(QASMstring)), QASMstring);
      
      % draw 
      [out] = H.draw(1, 'N');
      test.verifyEqual( out, 0 );
      [out] = H.draw(0, 'L');
      test.verifyTrue( isa(out, 'cell') );
      test.verifySize( out, [3, 1] );
      
      % TeX
      [out] = H.toTex(1, 'N');
      test.verifyEqual( out, 0 );
      [out] = H.toTex(0, 'L');
      test.verifyTrue( isa(out, 'cell') );
      test.verifySize( out, [1, 1] );

      % offset
      H.setOffset( 3 );
      test.verifyEqual( H.qubit, int64(3) );
      test.verifyEqual( H.offset, int64(3) );
      qubits = H.qubits;
      test.verifyEqual( qubits(1), int64(3) );
      [T,~] = evalc('H.toQASM(1)'); % capture output to std::out in T
      QASMstring = sprintf(['\n\nQASM output:\n\n' ...
                      'OPENQASM 2.0;\n' ...
                      'include "qelib1.inc";\n\n' ...
                      'qreg q[1];\n' ...
                      'x q[3];\n' ...
                    ]);
      test.verifyEqual(T(1:length(QASMstring)), QASMstring);
      
      % draw gate with offset
      [out] = H.draw(1, 'N');
      test.verifyEqual( out, 0 );
      [out] = H.draw(0, 'L');
      test.verifyTrue( isa(out, 'cell') );
      test.verifySize( out, [3, 1] );
      
      % TeX with offset
      [out] = H.toTex(1, 'N');
      test.verifyEqual( out, 0 );
      [out] = H.toTex(0, 'L');
      test.verifyTrue( isa(out, 'cell') );
      test.verifySize( out, [1, 1] );
      
      % handle
      HX = H.circuitHandle ;
      HX.setOffset( 3 );
      test.verifyEqual( H.qubit, int64(6) );
      CX = H.circuit ;
      CX.setOffset( 1 );
      test.verifyEqual( H.qubit, int64(6) );
      H.setCircuit(CX);
      test.verifyEqual( H.qubit, int64(4) );
      
      % operators == and ~=
      test.verifyTrue( H == Circ1 );
      test.verifyFalse( H ~= Circ1 );
      Z = qclab.qgates.PauliZ;
      Circ2 = qclab.QCircuit(1);
      Circ2.push_back(Z);
      test.verifyTrue( H ~= Circ2 );
      test.verifyFalse( H == Circ2 );
      X2 = qclab.qgates.PauliX;
      Circ3 = qclab.QCircuit(1);
      Circ3.push_back(X2);
      HX = qclab.HandleCircuit( Circ3 ) ;
      HZ = qclab.HandleCircuit( Circ2 ) ;
      test.verifyTrue( H == HX );
      test.verifyFalse( H ~= HX );
      test.verifyTrue( H ~= HZ );
      test.verifyFalse( H == HZ );
      
      % ctranspose
      Hp = H';
      Xp = X';
      test.verifyEqual(Hp.matrix, Xp.matrix );
    end

    function test_HandleCircuit_3qubits(test)
      X = qclab.qgates.PauliX ;
      Circ1 = qclab.QCircuit(3,1) ;
      Circ1.push_back(X) ;
      H = qclab.HandleCircuit(Circ1,2) ;

      test.verifyEqual( H.nbQubits, int64(3) );     % nbQubits
      test.verifyFalse( H.fixed );                   % fixed
      test.verifyFalse( H.controlled );             % controlled
      test.verifyEqual( H.qubit, int64(3) );        % qubit
      test.verifyEqual( H.offset, int64(2) );       % offset
      test.verifySameHandle( H.circuitHandle, Circ1 );     % circuitHandle

      % qubits
      qubits = H.qubits;
      test.verifyEqual( length(qubits), 3 );
      test.verifyEqual( qubits(1), int64(3) );
      test.verifyEqual( qubits(2), int64(4) );
      test.verifyEqual( qubits(3), int64(5) );
      
      % matrix
      test.verifyEqual( H.matrix, Circ1.matrix );
      
      % offset
       H.setOffset( 3 );
       test.verifyEqual( H.qubit, int64(4) );
       test.verifyEqual( H.offset, int64(3) );
       qubits = H.qubits;
       test.verifyEqual( qubits(1), int64(4) );
       test.verifyEqual( qubits(2), int64(5) );
       test.verifyEqual( qubits(3), int64(6) );
      
      % handle
       HX = H.circuitHandle ;
       HX.setOffset( 3 );
       test.verifyEqual( H.qubit, int64(6) );
       CX = H.circuit ;
       CX.setOffset( 1 );
       test.verifyEqual( H.qubit, int64(6) );
       H.setCircuit(CX);
       test.verifyEqual( H.qubit, int64(4) );
       
      % operators == and ~=
      test.verifyTrue( H == Circ1 );
      test.verifyFalse( H ~= Circ1 );
      Z = qclab.qgates.PauliZ;
      Circ2 = qclab.QCircuit(1);
      Circ2.push_back(Z);
      test.verifyTrue( H ~= Circ2 );
      test.verifyFalse( H == Circ2 );
      X2 = qclab.qgates.PauliX;
      Circ3 = qclab.QCircuit(3);
      Circ3.push_back(X2);
      HX = qclab.HandleCircuit( Circ3 ) ;
      HZ = qclab.HandleCircuit( Circ2 ) ;
      test.verifyTrue( H == HX );
      test.verifyFalse( H ~= HX );
      test.verifyTrue( H ~= HZ );
      test.verifyFalse( H == HZ );
      
      % ctranspose
       Hp = H';
       Circ1p = Circ1';
       test.verifyEqual(Hp.matrix, Circ1p.matrix );
    end
  end
end