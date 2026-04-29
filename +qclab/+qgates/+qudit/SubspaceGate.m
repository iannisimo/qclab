
classdef SubspaceGate < qclab.qgates.QGate1
  properties ( Access = protected )
    %> Child gate
    gate_ qclab.qgates.QGate1;
    %> Subspace
    subspace_ double;
  end
  
  methods

    function obj = SubspaceGate( gate, subspace, qubit )
      if nargin <= 1, subspace = [0,1]; end
      if nargin <= 2, qubit = 0; end
      assert(all(subspace >= 0));
      obj@qclab.qgates.QGate1(qubit);

      if isa(gate, 'qclab.qgates.QGate1')
        obj.gate_ = copy(gate);
        obj.gate_.setQubits(0);
      else
        obj.gate_ = gate( target );
        assert(isa(obj.gate_, 'qclab.qgates.QGate1'));
      end

      obj.subspace_ = subspace;
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      assert(d >= length(obj.subspace_));
      assert(d >= max(obj.subspace_) + 1);
      % TODO make issparse dependant on d and update this
      isSparse = false;
      submat = obj.gate_.matrix(length(obj.subspace_));
      mat = qclab.qId(1, isSparse, d);
      mat(obj.subspace_ + 1, obj.subspace_ + 1) = submat;
    end

    function [bool] = equals(obj, other)
      bool = isa(other, 'qclab.qgates.qudit.SubspaceGate') && ... % Same class
        obj.subspace_ == other.subspace_ && ... % Acting on the same subspace
        obj.gate_.equals(other.gate_); % Same gate
    end

    function [out] = toQASM(obj, fid, offset)
      assert(false, 'Unsupported');
    end


    function [out] = toTex(obj, fid, parameter, offset)
      assert(false, 'Unsupported');
    end

    


    

    
  end

  methods (Access = protected)
    function qd = isQudit(obj)
      qd = all(arrayfun(@(g) g.isQudit(), obj.gates_));
    end
  end

  methods (Static)
    function [bool] = fixed
      bool = true;
    end

    function [label] = label(obj, parameter, tex )
      label = 'X';
    end
  end

end
