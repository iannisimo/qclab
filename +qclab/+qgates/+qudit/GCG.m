%> GCG: Generator Controlled Gate
%> Wraps a single 1-qudit "generator" gate G into a controlled operation:
%> GCG|c,t> applies G^c to the target qudit t, i.e. the branch gate for
%> control value k is G^k. Generalizes constructions such as SUM, which
%> is recovered by taking G to be the cyclic shift (X) gate.
classdef GCG < qclab.qgates.QGate2
  properties ( Access = protected )
    %> Control qudit
    control_    int64
    %> Target qudit
    target_     int64
    %> Generator 1-qudit gate: branch k applies gate_^k to the target
    gate_       qclab.qgates.QGate1
  end

  methods

    function obj = GCG( control, target, gate )
      if nargin <= 0, control = 0; end
      if nargin <= 1, target = 1; end
      assert( qclab.isNonNegInteger(control) );
      assert( qclab.isNonNegInteger(target) );
      assert( control ~= target, 'control and target qudits cannot be the same' );
      assert( isa(gate, 'qclab.qgates.QGate1') );
      obj.control_ = control;
      obj.target_ = target;
      obj.gate_ = copy(gate);
      obj.gate_.setQubits(0);
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      isSparse = qclab.isSparse(obj.nbQubits, d);
      powers = obj.powers(d);
      if obj.control_ < obj.target_
        mat = blkdiag(powers{:});
      else
        mat = zeros(d^2, d^2);
        for k = 0:d-1
          mat = mat + kron(powers{k+1}, qclab.En(k, d, isSparse));
        end
      end
      if isSparse, mat = sparse(mat); end
    end

    function [current] = apply(obj, side, op, nbQubits, current, offset, d)
      if nargin <= 5, offset = 0; end
      if nargin <= 6, d = 2; end
      assert(nbQubits >= 2);
      isSparse = qclab.isSparse(nbQubits, d);
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
      if obj.control_ + 1 == obj.target_ || obj.target_ + 1 == obj.control_
        current = apply@qclab.qgates.QGate2( obj, side, op, nbQubits, ...
          current, offset, d);
        return
      end
      powers = obj.powers(d, op);
      s = abs(obj.target_ - obj.control_) + 1;
      Imid = qclab.qId(s-2, isSparse, d);
      mats = zeros(d^s, d^s);
      for k = 0:d-1
        Ek = qclab.En(k, d, isSparse);
        if obj.control_ < obj.target_
          mats = mats + kron(kron(Ek, Imid), powers{k+1});
        else
          mats = mats + kron(kron(powers{k+1}, Imid), Ek);
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
      if obj.control_ < obj.target_
        obj.control_ = qubits(1);
        obj.target_ = qubits(2);
      else
        obj.control_ = qubits(2);
        obj.target_ = qubits(1);
      end
    end

    function [control] = control(obj)
      control = obj.control_;
    end

    function [target] = target(obj)
      target = obj.target_;
    end

    %> Copy of the 1-qudit generator gate
    function [gate] = gate(obj)
      gate = copy(obj.gate_);
    end

    function [bool] = equals(obj, other)
      bool = isa(other, 'qclab.qgates.qudit.GCG') && ...
        (sign(obj.control_ - obj.target_) == sign(other.control_ - other.target_)) && ...
        obj.gate_.equals(other.gate_);
    end

    function [out] = toQASM(obj, fid, offset)
      assert(false, 'Unsupported');
    end

    function [out] = toTex(obj, fid, parameter, offset)
      assert(false, 'Unsupported');
    end

    function [label] = label(obj, parameter, tex)
      if nargin < 2, parameter = 'N'; end
      if nargin < 3, tex = false; end
      label = obj.gate_.label( parameter, tex );
    end

  end

  methods ( Access = protected )
    function qd = isQudit(obj)
      qd = obj.gate_.isQudit();
    end

    %> Cell array {G^0, G^1, ..., G^(d-1)} of powers of the generator
    %> gate's d x d matrix, optionally transposed/conjugated per `op`
    %> ('N', 'T' or 'C', matching apply's convention)
    function [powers] = powers(obj, d, op)
      if nargin <= 2, op = 'N'; end
      G = obj.gate_.matrix(d);
      if strcmp(op, 'T')
        G = G.';
      elseif ~strcmp(op, 'N')
        G = G';
      end
      powers = cell(1, d);
      powers{1} = eye(d);
      for k = 1:d-1
        powers{k+1} = powers{k} * G;
      end
    end

    %> Override copyElement function to allow for correct deep copy of
    %> handle property
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable( obj );
      cp.gate_ = obj.gate();
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
