%> MVCG: Multi Valued Multi-Control Gate
classdef MVMCG < qclab.QObject
  properties ( Access = protected )
    %> Control qudit
    controls_    int64
    %> Target qudit
    target_     int64
    %> Number of qudits for the gate
    nbQubits_   int64
    %> Ordered list of gates to apply at each control value
    gates_(1,:) qclab.qgates.QGate1
  end
  
  methods

    function obj = MVMCG( controls, target, gates )
      assert(length(controls) == length(unique(controls)))
      assert(all(controls >= 0));
      assert(target >= 0);
      assert(control ~= target, 'control and target qudits cannot be the same');
      assert(isa(gates, 'qclab.qgates.QGate1'));
      [obj.controls_, sortIdx] = sort(controls);
      obj.target_ = target;
      % TODO handle @gates
      obj.gates_ = copy(gates(sortIdx));
      arrayfun(@(g) g.setQubits(0), obj.gates_);
      obj.nbQubits_ = length(controls) + 1;
    end

    function [mat] = matrix(obj, d)
      if nargin <= 1, d = 2; end
      assert(length(obj.gates_) == d^length(obj.controls_));
      % TODO make issparse dependant on d and update this
      isSparse = false;
      if obj.controls_ < obj.target_
        matrices = arrayfun(@(g) g.matrix(d), obj.gates_, 'UniformOutput', false);
        mat = blkdiag(matrices{:});
      else
        mat = zeros(d^obj.nbQubits_, d^obj.nbQubits_);
        for i = 1:length(obj.gates_)
          mat = mat + kron(obj.gates_(i).matrix(d), qclab.En(i-1, d^length(obj.controls_), isSparse));
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
      if (obj.controls_(end) + 1 == obj.target_) || (obj.target_ + 1 == obj.controls_(1))
        % operation
        if strcmp(op, 'N') % normal
          mat = obj.matrix(d);
        elseif strcmp(op, 'T') % transpose
          mat = obj.matrix(d).';
        else % conjugate transpose
          mat = obj.matrix(d)';
        end
        % kron( Ileft, mat2, Iright)
        if (nbQubits == obj.nbQubits_)
          matn = mat ;
        elseif ( qubits(1) == 0 )
          matn = kron(mat, qclab.qId(nbQubits-obj.nbQubits_, isSparse, d)) ;
        elseif ( qubits(end) == nbQubits-1)
          matn = kron(qclab.qId(nbQubits-obj.nbQubits_, isSparse, d), mat);
        else
          matn = kron(kron(qclab.qId(qubits(1), isSparse, d), mat), ...
            qclab.qId(nbQubits-qubits(end)-1, isSparse, d)) ;
        end
        current = qclab.applyGateTo( current, matn, side ) ;
        return
      end
      s = qubits(end) - qubits(1) + 1;
      mats = zeros(d^s, d^s);
      for i = 1:length(obj.gates_)
        C = repmat(-1, s);
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
      qubits = sort([obj.controls_, obj.target_]);
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
