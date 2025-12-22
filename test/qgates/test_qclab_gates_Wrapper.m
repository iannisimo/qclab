classdef test_qclab_gates_Wrapper < matlab.unittest.TestCase
  methods (Test)

    function test_Hadamard(testCase)
      circuit = qclab.QCircuit(1);
      circuit.H(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.Hadamard(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_CNOT(testCase)
      circuit = qclab.QCircuit(2);
      circuit.CNOT(0, 1);
      gate = circuit.objects(1);
      expected = qclab.qgates.CNOT(0, 1);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_CPhase(testCase)
      theta = pi;
      circuit = qclab.QCircuit(2);
      circuit.CPhase(0, 1, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.CPhase(0, 1, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_CRX(testCase)
      theta = pi/2;
      circuit = qclab.QCircuit(2);
      circuit.CRX(0, 1, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.CRotationX(0, 1, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_CRY(testCase)
      theta = pi/3;
      circuit = qclab.QCircuit(2);
      circuit.CRY(0, 1, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.CRotationY(0, 1, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_CRZ(testCase)
      theta = pi/4;
      circuit = qclab.QCircuit(2);
      circuit.CRZ(0, 1, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.CRotationZ(0, 1, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_CU2(testCase)
      phi = pi/5; lambda = pi/6;
      circuit = qclab.QCircuit(2);
      circuit.CU2(0, 1, 1, phi, lambda);
      gate = circuit.objects(1);
      expected = qclab.qgates.CU2(0, 1, 1, phi, lambda);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.phi, phi, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.lambda, lambda, 'AbsTol', 1e-14);
    end

    function test_CU3(testCase)
      theta = 1; phi = pi/5; lambda = pi/6;
      circuit = qclab.QCircuit(2);
      circuit.CU3(0, 1,1, theta, phi, lambda);
      gate = circuit.objects(1);
      expected = qclab.qgates.CU3(0, 1,1, theta, phi, lambda);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.phi, phi, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.lambda, lambda, 'AbsTol', 1e-14);
    end

    function test_CY(testCase)
      circuit = qclab.QCircuit(2);
      circuit.CY(0, 1);
      gate = circuit.objects(1);
      expected = qclab.qgates.CY(0, 1);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_CZ(testCase)
      circuit = qclab.QCircuit(2);
      circuit.CZ(0, 1);
      gate = circuit.objects(1);
      expected = qclab.qgates.CZ(0, 1);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_iSWAP(testCase)
      circuit = qclab.QCircuit(2);
      circuit.iSWAP(0, 1);
      gate = circuit.objects(1);
      expected = qclab.qgates.iSWAP(0, 1);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_MG(testCase)
      circuit = qclab.QCircuit(2);
      U = [0 1; 1 0];
      circuit.MG([0], U);
      gate = circuit.objects(1);
      expected = qclab.qgates.MatrixGate([0], U);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.matrix, U, 'AbsTol', 1e-14);
    end

    function test_MCMG(testCase)
      circuit = qclab.QCircuit(3);
      U = eye(2);
      circuit.MCMG([0 1], 2, U);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCMatrixGate([0 1], 2, U);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_MCRX(testCase)
      theta = pi/2;
      circuit = qclab.QCircuit(3);
      circuit.MCRX([0 1], 2, [1 0], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCRotationX([0 1], 2, [1 0], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_MCRY(testCase)
      theta = pi/3;
      circuit = qclab.QCircuit(3);
      circuit.MCRY([0 1], 2, [1 1], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCRotationY([0 1], 2, [1 1], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_MCRZ(testCase)
      theta = pi/4;
      circuit = qclab.QCircuit(3);
      circuit.MCRZ([0 1], 2, [0 1], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCRotationZ([0 1], 2, [0 1], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_MCX(testCase)
      circuit = qclab.QCircuit(3);
      circuit.MCX([0 1], 2);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCX([0 1], 2);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_MCY(testCase)
      circuit = qclab.QCircuit(3);
      circuit.MCY([0 1], 2);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCY([0 1], 2);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_MCZ(testCase)
      circuit = qclab.QCircuit(3);
      circuit.MCZ([0 1], 2);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCZ([0 1], 2);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_PauliX(testCase)
      circuit = qclab.QCircuit(1);
      circuit.X(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.PauliX(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_PauliY(testCase)
      circuit = qclab.QCircuit(1);
      circuit.Y(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.PauliY(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_PauliZ(testCase)
      circuit = qclab.QCircuit(1);
      circuit.Z(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.PauliZ(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_Phase(testCase)
      theta = pi/8;
      circuit = qclab.QCircuit(1);
      circuit.Phase(0, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase(0, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_Phase45(testCase)
      circuit = qclab.QCircuit(1);
      circuit.Phase45(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase45(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_T(testCase)
      circuit = qclab.QCircuit(1);
      circuit.T(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase45(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_Phase90(testCase)
      circuit = qclab.QCircuit(1);
      circuit.Phase90(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase90(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_S(testCase)
      circuit = qclab.QCircuit(1);
      circuit.S(0);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase90(0);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_RX(testCase)
      theta = pi/2;
      circuit = qclab.QCircuit(1);
      circuit.RX(0, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationX(0, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_RXX(testCase)
      theta = pi/2;
      circuit = qclab.QCircuit(2);
      circuit.RXX([0 1], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationXX([0 1], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_RY(testCase)
      theta = pi/3;
      circuit = qclab.QCircuit(1);
      circuit.RY(0, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationY(0, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_RYY(testCase)
      theta = pi/3;
      circuit = qclab.QCircuit(2);
      circuit.RYY([0 1], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationYY([0 1], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_RZ(testCase)
      theta = pi/4;
      circuit = qclab.QCircuit(1);
      circuit.RZ(0, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationZ(0, theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_RZZ(testCase)
      theta = pi/4;
      circuit = qclab.QCircuit(2);
      circuit.RZZ([0 1], theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationZZ([0 1], theta);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
    end

    function test_SWAP(testCase)
      circuit = qclab.QCircuit(2);
      circuit.SWAP(0, 1);
      gate = circuit.objects(1);
      expected = qclab.qgates.SWAP(0, 1);
      testCase.verifyTrue(gate.equals(expected));
    end

    function test_U2(testCase)
      phi = pi/3; lambda = pi/4;
      circuit = qclab.QCircuit(1);
      circuit.U2(0, phi, lambda);
      gate = circuit.objects(1);
      expected = qclab.qgates.U2(0, phi, lambda);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.phi, phi, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.lambda, lambda, 'AbsTol', 1e-14);
    end

    function test_U3(testCase)
      theta = pi/5; phi = pi/6; lambda = pi/7;
      circuit = qclab.QCircuit(1);
      circuit.U3(0, theta, phi, lambda);
      gate = circuit.objects(1);
      expected = qclab.qgates.U3(0, theta, phi, lambda);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, theta, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.phi, phi, 'AbsTol', 1e-14);
      testCase.verifyEqual(gate.lambda, lambda, 'AbsTol', 1e-14);
    end

    function test_CNOT_With_ControlState0(testCase)
      circuit = qclab.QCircuit(2);
      circuit.CNOT(0, 1, 0);
      gate = circuit.objects(1);
      expected = qclab.qgates.CNOT(0, 1, 0);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.controlState, int64(0));
    end

    function test_Multiple_Gates_In_Circuit(testCase)
      circuit = qclab.QCircuit(2);
      circuit.H(0);
      circuit.CNOT(0, 1);
      objects = circuit.objects;
      testCase.verifyEqual(numel(objects), 2);
      testCase.verifyTrue(objects(1).equals(qclab.qgates.Hadamard(0)));
      testCase.verifyTrue(objects(2).equals(qclab.qgates.CNOT(0, 1)));
    end

    function test_Empty_MatrixGate_Label(testCase)
      circuit = qclab.QCircuit(1);
      U = [0 1; 1 0];
      circuit.MG(0, U);
      gate = circuit.objects(1);
      expected = qclab.qgates.MatrixGate(0, U);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.label, ' U ');
    end

    function test_Custom_Label_MatrixGate(testCase)
      circuit = qclab.QCircuit(1);
      U = [1 0; 0 -1];
      circuit.MG(0, U, 'Z');
      gate = circuit.objects(1);
      expected = qclab.qgates.MatrixGate(0, U, 'Z');
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.label, ' Z ');
    end

    function test_MCMG_With_ControlStates(testCase)
      circuit = qclab.QCircuit(3);
      U = [1 0; 0 -1];
      circuit.MCMG([0 1], 2, U, [1 0]);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCMatrixGate([0 1], 2, U, [1 0]);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.controlStates, int64([1 0]));
    end

    function test_RX_Default_Qubit(testCase)
      circuit = qclab.QCircuit(1);
      circuit.RX();
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationX(0, 0);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.qubits, int64(0));
    end

    function test_RY_With_QAngle(testCase)
      circuit = qclab.QCircuit(1);
      angle = qclab.QAngle(pi/2);
      circuit.RY(0, angle);
      gate = circuit.objects(1);
      expected = qclab.qgates.RotationY(0, pi);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, pi, 'AbsTol', 1e-14);
    end


    function test_U2_With_Default_Args(testCase)
      circuit = qclab.QCircuit(1);
      circuit.U2();
      gate = circuit.objects(1);
      expected = qclab.qgates.U2(0, 0, 0);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.qubits, int64(0));
    end

    function test_MCZ_With_ControlState0(testCase)
      circuit = qclab.QCircuit(3);
      circuit.MCZ([0 1], 2, [0 1]);
      gate = circuit.objects(1);
      expected = qclab.qgates.MCZ([0 1], 2, [0 1]);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.controlStates, int64([0 1]));
    end

    function test_Phase_With_Angle_Object(testCase)
      circuit = qclab.QCircuit(1);
      theta = qclab.QAngle(pi/8);
      circuit.Phase(0, theta);
      gate = circuit.objects(1);
      expected = qclab.qgates.Phase(0, pi/8);
      testCase.verifyTrue(gate.equals(expected));
      testCase.verifyEqual(gate.theta, pi/8, 'AbsTol', 1e-14);
    end
  
  end
end