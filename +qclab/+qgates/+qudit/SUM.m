%> SUM: qudit SUM gate
%> Generalization of CNOT to qudits: SUM|c,t> = |c, (c+t) mod d>
classdef SUM < qclab.qgates.QGate2
  properties ( Access = protected )
    %> Control qudit
    control_    int64
    %> Target qudit
    target_     int64
  end

  methods

    function obj = SUM( control, target )
      if nargin <= 0, control = 0; end
      if nargin <= 1, target = 1; end
      assert( qclab.isNonNegInteger(control) );
      assert( qclab.isNonNegInteger(target) );
      assert( control ~= target, 'control and target qudits cannot be the same' );
      obj.control_ = control;
      obj.target_ = target;
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      % TODO make issparse dependent on d and update this
      isSparse = false;
      Id = eye(d);
      blocks = arrayfun(@(k) circshift(Id, k), 0:d-1, 'UniformOutput', false);
      if obj.control_ < obj.target_
        mat = blkdiag(blocks{:});
      else
        mat = zeros(d^2, d^2);
        for k = 0:d-1
          mat = mat + kron(blocks{k+1}, qclab.En(k, d, isSparse));
        end
      end
      if isSparse, mat = sparse(mat); end
    end

    function [current] = apply(obj, side, op, nbQubits, current, offset, d)
      if nargin <= 5, offset = 0; end
      if nargin <= 6, d = 2; end
      isSparse = qclab.isSparse(nbQubits);
      assert( nbQubits >= 2 );
      if isa(current, 'double')
        if strcmp(side,'L') % left
          assert( size(current,2) == d^nbQubits);
        else % right
          assert( size(current,1) == d^nbQubits);
        end
      else
        assert( length(current.states{1}) == d^nbQubits )
      end
      qubits = obj.qubits + offset;
      assert( qubits(1) < nbQubits ); assert( qubits(2) < nbQubits );
      % nearest neighbor qudits
      if qubits(1) + 1 == qubits(2)
        current = apply@qclab.qgates.QGate2( obj, side, op, nbQubits, ...
          current, offset, d);
        return
      end
      Id = eye(d);
      s = qubits(2) - qubits(1) + 1;
      Imid = qclab.qId(s-2, isSparse, d);
      mats = zeros(d^s, d^s);
      for k = 0:d-1
        block = circshift(Id, k);
        if strcmp(op, 'T') % transpose
          block = block.';
        elseif ~strcmp(op, 'N') % conjugate transpose
          block = block';
        end
        Ek = qclab.En(k, d, isSparse);
        if obj.control_ < obj.target_
          mats = mats + kron(kron(Ek, Imid), block);
        else
          mats = mats + kron(kron(block, Imid), Ek);
        end
      end
      if isSparse, mats = sparse(mats); end
      if ( qubits(1) == 0 && qubits(2) == nbQubits - 1)
        matn = mats;
      elseif ( qubits(1) == 0 )
        matn = kron(mats, qclab.qId(nbQubits - s, isSparse, d));
      elseif ( qubits(2) == nbQubits - 1 )
        matn = kron(qclab.qId(nbQubits - s, isSparse, d), mats);
      else
        matn = kron(kron(qclab.qId(qubits(1), isSparse, d), mats), ...
          qclab.qId(nbQubits - qubits(2) - 1, isSparse, d));
      end
      current = qclab.applyGateTo(current, matn, side );
    end

    function [qubit] = qubit(obj)
      qubit = min(obj.control_, obj.target_);
    end

    function [qubits] = qubits(obj)
      qubits = sort([obj.control_, obj.target_]);
    end

    function setQubits(obj, qubits)
      assert( qclab.isNonNegIntegerArray(qubits) );
      assert( qubits(1) ~= qubits(2) );
      obj.control_ = qubits(1);
      obj.target_ = qubits(2);
    end

    function [control] = control(obj)
      control = obj.control_;
    end

    function [target] = target(obj)
      target = obj.target_;
    end

    function [bool] = equals(obj, other)
      bool = isa(other, 'qclab.qgates.qudit.SUM') && ...
        (sign(obj.control_ - obj.target_) == sign(other.control_ - other.target_));
    end

    function [out] = toQASM(obj, fid, offset)
      assert(false, 'Unsupported');
    end

    function [out] = toTex(obj, fid, parameter, offset)
      assert(false, 'Unsupported');
    end

    function [label] = label(obj, parameter, tex )
      label = char(8853); % oplus symbol
    end

  end

  methods (Access = protected)
    function qd = isQudit(~)
      qd = true;
    end
  end

  methods (Static)
    function [bool] = controlled
      bool = false;
    end

    function [bool] = fixed
      bool = true;
    end
  end

end
