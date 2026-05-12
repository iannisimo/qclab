classdef ControlledGate < qclab.qgates.QControlledGate2
  
  properties (Access = protected)
    %> Property storing the 1-qubit gate of this 2-qubit controlled gate.
    gate_ qclab.QObject
  end
  
  methods
    
    function obj = ControlledGate(gate, control, target, controlState)
      % set defaults
      if nargin <= 1, control = 0; end
      if nargin <= 2, target = 1; end
      if nargin <= 3, controlState = 1; end
      obj@qclab.qgates.QControlledGate2(control, controlState );

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

      assert( qclab.isNonNegInteger(target) ); 
      assert(control ~= obj.target) ;
      
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      assert(false, 'TODO: unimplemented')
      if nargin == 2, offset = 0; end
    end
    
    % equals
    function [bool] = equals(obj, other)
      bool = false;
      if isa(other,'qclab.qgates.ControlledGate') && ...
          (obj.controlState == other.controlState)
        bool = ((obj.control < obj.target) && (other.control < other.target))...
            || ((obj.control > obj.target) && (other.control > other.target)) && ...
          obj.gate_.equals(other.gate_);
      end
    end
    
    % target
    function [target] = target(obj)
      target = obj.gate_.qubit ;
    end
    
    %> Copy of 1-qubit gate of controlled-gate
    function [gate] = gate(obj)
      gate = copy(obj.gate_);
    end
    
      
    % setTarget
    function setTarget(obj, target)
      assert( qclab.isNonNegInteger(target) ) ; 
      assert( target ~= obj.control() ) ;
      obj.gate_.setQubit( target );
    end
    %
    % label for draw and tex function
    function [label] = label(obj, parameter, tex )
      if nargin < 2, parameter = 'N'; end
      if nargin < 3, tex = false; end
      label = obj.gate_.label( parameter, tex );
    end
    
  end
  
  methods (Static)
    
    function [bool] = fixed
      bool = obj.gate_.fixed;
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
