% PauliX - 1-qubit Pauli-X gate for quantum circuits
% The PauliX class implements a 1-qubit Pauli-X gate. The Pauli-X gate 
% performs a bit flip, swapping the |0⟩ and |1⟩ states. It is the quantum 
% analog of the classical NOT gate.
%
% The matrix representation of the Pauli-X gate is:
%   X = [0  1; 
%        1  0]
%
% Creation
%   Syntax
%     X = qclab.qgates.PauliX(qubit)
%
%   Input Arguments
%     qubit - qubit to which the Pauli-X gate is applied
%             non-negative integer, (default: 0)
%
%   Output:
%     X - A quantum object of type `PauliX`, representing the 1-qubit
%         Pauli-X gate on qubit `qubit`.
%
% Example:
%   Create a Pauli-X gate object acting on qubit 0:
%     X = qclab.qgates.PauliX(0);

%> @file PauliX.m
%> @brief Implements Pauli-X class.
% ==============================================================================
%> @class PauliX
%> @brief 1-qubit Pauli-X gate.
%>
%> 1-qubit Pauli-X gate with matrix representation:
%>
%> \f[\begin{bmatrix} 0 & 1\\ 
%>                    1 & 0 \end{bmatrix}\f]
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef PauliX < qclab.qgates.QGate1
  properties (Access = protected)
    %> true if this gate is the adjoint (differs from the gate for d > 2)
    adjoint_(1,1) logical = false
  end

  methods (Static)
    % fixed
    function [bool] = fixed
      bool = true;
    end
  end
  
  methods
    % matrix
    function [mat] = matrix(obj, d)
      if nargin < 2
        d = 2;
      end
      if d == 2
        mat = [0 1; 1 0];
        return
      end
      mat = circshift(eye(d), 1, 1);
      if obj.adjoint_, mat = mat'; end
    end
    
    % label for draw and tex function
    function [label] = label(obj, parameter, tex )
      label = 'X';
      if obj.adjoint_, label = 'X'''; end
    end

    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end
      qclab.IO.qasmPauliX( fid, obj.qubit + offset );
      out = 0;
    end
    
    % equals
    function [bool] = equals(obj,other)
      bool = isa(other, 'qclab.qgates.PauliX') && ...
        obj.adjoint_ == other.adjoint_;
    end

    % ctranspose: for d > 2 the gate is not self-adjoint
    function objprime = ctranspose( obj )
      objprime = copy( obj );
      objprime.adjoint_ = ~obj.adjoint_;
    end
  end

  methods ( Access = protected )
    function qd = isQudit(~)
      qd = true;
    end
  end
end % PauliX
