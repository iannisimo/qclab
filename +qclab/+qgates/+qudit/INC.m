classdef INC < qclab.qgates.QGate1 & qclab.QAdjustable
  properties ( Access = protected )
    sumval_(1,1) int32
  end
  
  methods

    function obj = INC( qubit, sumval, fixed )
      if nargin <= 0, qubit = 0; end
      if nargin <= 1, sumval = 1; end
      if nargin <= 2, fixed = false; end
      obj@qclab.qgates.QGate1(qubit);
      obj@qclab.QAdjustable(fixed); 
      obj.sumval_ = sumval;
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      % TODO make issparse dependent on d and update this
      isSparse = true;
      mat = qclab.qId(1, isSparse, d);
      mat = circshift(mat, obj.sumval_);
    end

    function update(obj, sumval)
      assert( ~obj.fixed );
      obj.sumval_ = sumval;
    end

    function [bool] = equals(obj, other)
      bool = isa(other, 'qclab.qgates.qudit.INC') && (obj.sumval_ == other.sumval_);
    end

    function [out] = toQASM(obj, fid, offset)
      assert(false, 'Unsupported');
    end


    function [out] = toTex(obj, fid, parameter, offset)
      assert(false, 'Unsupported');
    end

    function [label] = label(obj, parameter, tex )
      label = [0x2295, int2str(obj.sumval_)];
    end

  end

  methods (Access = protected)
    function qd = isQudit(~)
      qd = true;
    end
  end

end
