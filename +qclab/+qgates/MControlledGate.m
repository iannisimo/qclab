classdef MControlledGate < qclab.qgates.QMultiControlledGate
  properties (Access = protected)
    %> Property storing the 1-qubit gate of this multi-qubit controlled 
    %> gate.
    gate_ qclab.QObject
  end
  
  methods
    
    function obj = MControlledGate(gate, controls, target, controlStates)
      if nargin <= 3, controlStates = ones(size(controls)) ; end
      obj@qclab.qgates.QMultiControlledGate(controls, controlStates );

      if isa(gate, 'qclab.qgates.QGate1')
        obj.gate_ = copy(gate);
        obj.gate_.setQubits(target);
      elseif isa(gate, 'qclab.QObject')
        assert(gate.nbQubits == 1);
        obj.gate_ = copy(gate);
        obj.gate_.setQubits(target);
      else
        obj.gate_ = gate( target );
        assert(isa(obj.gate_, 'qclab.qgates.QGate1'));
      end

      assert( qclab.isNonNegInteger(target) ) ; 
      assert( ~(any(controls(:) == target) ) ) ;
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      assert(false, 'TODO: unimplemented')
    end
    
    % equals
    function [bool] = equals(obj, other)
      bool = false ;
      if isa(other,'qclab.qgates.MControlledGate') && ...
         all(obj.controlStates == other.controlStates)
        bool = (sum(obj.controls < obj.targets) == ...
                sum(other.controls < other.targets)) && ...
               (sum(obj.controls > obj.targets) == ...
                sum(other.controls > other.targets)) && ...
                obj.gate_.equals(other.gate_) ;
      end
    end
    
     % targets
    function [targets] = targets(obj)
      targets = obj.gate_.qubit ;
    end
    
    % setTargets
    function setTargets(obj, target)
      assert( qclab.isNonNegInteger(target) ) ; 
      controls = obj.controls() ;
      assert( ~(any(controls(:) == target) ) ) ;
      obj.gate_.setQubit( target );
    end

    %> Copy of 1-qubit gate of multi-controlled gate
    function [gate] = gate(obj)
      gate = copy(obj.gate_);
    end
    
    % label for draw and tex function
    function [label] = label(obj, parameter, tex )
      if nargin < 2, parameter = 'N'; end
      if nargin < 3, tex = false; end
      label = obj.gate_.label( parameter, tex );
    end
    
    % setQubits
    function setQubits(obj, qubits)
      assert( qclab.isNonNegIntegerArray(qubits) ) ;
      assert( length(qubits) == length(unique(qubits)) );
      obj.controls_ = sort(qubits(1:end-1)) ;
      obj.setTargets( qubits(end) ) ;
    end
  end
  
  methods (Static)
    
    function [bool] = fixed
      bool = true;
    end
    
  end
  
  methods ( Access = protected )
    
    %> @brief Override copyElement function to allow for correct deep copy of
    %> handle property
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable( obj );
      cp.gate_ = obj.gate() ;
    end
    
  end
  
end
