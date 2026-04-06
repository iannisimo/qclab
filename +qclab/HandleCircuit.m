% HandleCircuit - Handle for quantum circuits
%
% HandleCircuit stores a handle to a quantum circuit together with a qubit
% offset. This allows the same circuit to be reused and applied to different
% qubits by shifting its qubit indices.
%
% If the referenced circuit itself already uses a qubit offset, the effective
% qubit indices are obtained by adding both offsets.
%
% Creation
%   Syntax
%     H = qclab.HandleCircuit(circuit)
%     H = qclab.HandleCircuit(circuit, offset)
%
%   Input Arguments
%     circuit - quantum circuit handle
%               object of type qclab.QCircuit
%
%     offset  - qubit offset applied to the circuit
%               non-negative integer (default: 0)
%
%   Output
%     H - quantum object of type HandleCircuit, representing a handle to
%         the given circuit with the specified qubit offset
%
% Example
%   Create a handle circuit and apply it to different qubit positions:
%     C = qclab.QCircuit(2);
%     H = qclab.HandleCircuit(C, 1);
%     H = qclab.HandleCircuit(C, 2);

%> @file HandleCircuit.m
%> @brief Implements HandleCircuit class
% ==============================================================================
%> @class HandleCircuit
%> @brief Quantum object storing a handle to another QCircuit.
%
%> 
% (C) Copyright Sophia Keip 2025.  
% ==============================================================================
classdef HandleCircuit < qclab.QObject & qclab.QAdjustable
  properties (Access = protected)
    %> Qubit offset of this handle circuit.
    offset_(1,1) int64
    %> Gate handle of this handle circuit.
    circuit_(1,1) %qclab.QCircuit
  end
  
  methods
    % Class constructor  =======================================================
    %> @brief Constructor for HandleCircuit objects
    %>
    %> Constructs a handle circuit from the given circuit handle `circuit`
    %> and qubit offset `offset`. The default value of `offset` is 0.
    % ==========================================================================
    function [obj] = HandleCircuit( circuit, offset )
      if nargin == 1, offset = 0; end
      obj.offset_ = offset;
      obj.circuit_ = circuit;
    end

    %> @brief Returns the qubit offset of this handle circuit.
    function [offset] = offset(obj)
      offset = obj.offset_ ;
    end

    %> @brief Sets the qubit offset of this handle circuit.
    function setOffset(obj, offset)
      obj.offset_ = offset ;
    end
    
    %> @brief Returns a handle to the circuit of this handle circuit.
    function [circuit] = circuitHandle(obj)
      circuit = obj.circuit_ ;
    end
    
    %> @brief Returns a copy of the circuit of this handle circuit.
    function [circuit] = circuit(obj)
      circuit = copy(obj.circuit_) ;
    end
    
    %> @brief Sets `circuit` as the new handle
    function setCircuit(obj, circuit)
      obj.circuit_ = circuit;
    end
    
    % nbQubits
    function [nbQubits] = nbQubits(obj)
      nbQubits = obj.circuit_.nbQubits ;
    end

    % fixed
    function [bool] = fixed(obj)
      bool = obj.circuit_.fixed ;
    end

    % controlled
    function [bool] = controlled(obj)
      bool = obj.circuit_.controlled ;
    end

    % qubit
    function [qubit] = qubit(obj)
      qubit = obj.circuit_.qubit + obj.offset_ ;
    end

    % setQubit
    function setQubit(~)
      assert( false );
    end
    
    % qubits
    function [qubits] = qubits(obj)
      qubits = obj.circuit_.qubits + obj.offset_ ;
    end 

    % setQubits
    function setQubits(~)
      assert( false );
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end        
      out = obj.circuit_.toQASM( fid, obj.offset_ + offset );
    end

    % equals
    function [bool] = equals(obj, other)
      if isa(other, 'qclab.HandleCircuit')
        bool = (other.circuit_ == obj.circuit_);
      else
        bool = (other == obj.circuit_ );
      end
    end

    % draw
    function [out] = draw(obj, fid, parameter, offset)
      if nargin < 2, fid = 1; end
      if nargin < 3, parameter = 'N'; end
      if nargin < 4, offset = 0; end
      out = obj.circuit_.draw(fid, parameter, obj.offset_ + offset );
    end

    % toTex
    function [out] = toTex(obj, fid, parameter, offset)
      if nargin < 2, fid = 1; end
      if nargin < 3, parameter = 'N'; end
      if nargin < 4, offset = 0; end
      out = obj.circuit_.toTex(fid, parameter, obj.offset_ + offset );
    end

    % apply
    function [current] = apply(obj, side, op, nbQubits, current, offset)
      if nargin == 5, offset = 0; end
      current = obj.circuit_.apply( side, op, nbQubits, current, ...
                                    offset + obj.offset_ ) ;
    end

    % matrix 
    function [mat] = matrix(obj)
      mat = obj.circuit_.matrix ;
    end

    % ctranspose
    function [objprime] = ctranspose( obj )
      objprime = copy( obj );
      objprime.setCircuit( ctranspose( obj.circuit_ ) );
    end
  end
  
  methods ( Access = protected )
    
    %> @brief Override copyElement function to allow for correct deep copy of
    %> handle property.
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable( obj );
      cp.circuit_ = obj.circuit() ;
    end

  end

end %HandleCircuit