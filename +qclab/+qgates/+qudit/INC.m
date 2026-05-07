classdef INC < qclab.qgates.QGate1 & qclab.QAdjustable
  properties ( Access = protected )
    sumval_(1,1) int32
  end
  
  methods

    function obj = INC( qubit, sumval )
      if nargin <= 1, qubit = 0; end
      if nargin <= 2, val = 1; end
      obj@qclab.qgates.QGate1(qubit);
      obj.sumval_ = sumval;
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      % TODO make issparse dependent on d and update this
      isSparse = false;
      mat = qclab.qId(1, isSparse, d);
      mat = circshift(mat, obj.sumval_);
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

  methods (Static)
    function [bool] = fixed
      bool = true;
    end
  end

end
