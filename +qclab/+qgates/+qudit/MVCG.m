%> MVCG: Multi Valued Control Gate
%> sect. 2.2.3 https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.589504/full
classdef MVCG < qclab.qgates.QGate2
  properties ( Access = protected )
    %> Control qudit
    control_    int64
    %> Target qudit
    target_     int64
    %> Ordered list of gates to apply at each control value
    gates_(1,:) qclab.qgates.QGate1
  end
  
  methods

    function obj = MVCG( control, target, gates )
      assert(control >= 0);
      assert(target >= 0);
      assert(control ~= target, 'control and target qudits cannot be the same');
      assert(isa(gates, 'qclab.qgates.QGate1'));
      obj.control_ = control;
      obj.target_ = target;
      obj.gates_ = copy(gates);
      arrayfun(@(g) g.setQubits(0), obj.gates_);
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      assert(length(obj.gates_) == d);
      % TODO make issparse dependant on d and update this
      isSparse = false;
      if obj.control_ < obj.target_
        matrices = arrayfun(@(g) g.matrix(d), obj.gates_, 'UniformOutput', false);
        mat = blkdiag(matrices{:});
      else
        mat = zeros(d^2, d^2);
        for i = 1:length(obj.gates_)
          mat = mat + kron(obj.gates_(i).matrix(d), qclab.En(i-1, d, isSparse));
        end
      end
      if isSparse, mat = sparse(mat); end
    end

    function [current] = apply(obj, side, op, nbQubits, current, offset, d)
      if nargin <= 5, offset = 0; end
      if nargin <= 6, d = 2; end
      assert(nbQubits >= 2);
      isSparse = qclab.isSparse(nbQubits);
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
      % nearest neighbor qubits
      if obj.control_ + 1 == obj.target_
        current = apply@qclab.qgates.QGate2( obj, side, op, nbQubits, ...
          current, offset, d);
        return
      end
      s = obj.target_ - obj.control_ + 1;
      Imid = qclab.qId(s-2, isSparse, d);
      mats = zeros(d^s, d^s);
      for i = 1:length(obj.gates_)
        if obj.control_ < obj.target_
          mats = mats + kron(kron(qclab.En(i-1, d, isSparse), Imid), obj.gates_(i).matrix(d));
        else
          mats = mats + kron(kron(obj.gates_(i).matrix(d), Imid), qclab.En(i-1, d, isSparse));
        end
      end
      if ( obj.control_ == 0 && obj.target_ == nbQubits - 1)
        matn = mats;
      elseif ( obj.control_ == 0 )
        matn = kron(mats, qclab.qId(nbQubits - s, isSparse, d));
      elseif ( obj.target_ == nbQubits - 1 )
        matn = kron(qclab.qId(nbQubits - s, isSparse, d), mats);
      else
        matn = kron(kron(qclab.qId(obj.control_, isSparse, d),mats),...
          qclab.qId(nbQubits - obj.target_ - 1, isSparse, d));
      end
      % apply
      current = qclab.applyGateTo(current, matn, side ) ;
    end

    function [bool] = equals(obj, other)
      bool = isa(other, 'qclab.qgates.qudit.MVCG') && ... % Same class
        (sign(obj.control_ - obj.target_) == sign(other.control_ - other.target_)) && ... % same direction
        length(obj.gates_) == length(other.gates_); % same number of gates / controls
      if ~bool, return; end
      for i = 1:length(obj.gates_)
        if ~obj.gates_(i).equals(other.gates_(i)), bool = false; return; end % same gates
      end
    end

    function [out] = toQASM(obj, fid, offset)
      assert(false, 'Unsupported');
    end


    function [out] = toTex(obj, fid, parameter, offset)
      assert(false, 'Unsupported');
    end

    function [qubit] = qubit(obj)
      qubit = min(obj.control_, obj.target_) ;
    end


    function [qubits] = qubits(obj)
      qubits = sort([obj.control_, obj.target_]);
    end

    function setQubits(obj, qubits)
      assert( qclab.isNonNegIntegerArray(qubits) ) ;
      assert( qubits(1) ~= qubits(2) ) ;
      obj.control_ = qubits(1) ;
      obj.target_ = qubits(2) ;
    end
  end

  methods (Access = protected)
    function qd = isQudit(obj)
      qd = all(arrayfun(@(g) g.isQudit(), obj.gates_));
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
