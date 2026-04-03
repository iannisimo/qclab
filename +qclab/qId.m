%> @file qId.m
%> @brief Implements quantum Identity matrix
% ==============================================================================
%> @brief Identity matrix on n qubits
%
%> @param n number of qubits
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.
% ==============================================================================
function [I] = qId(n,issparse,d)
if nargin < 2; issparse = false; end
if nargin < 3; d = 2; end
if issparse
  I = speye(d^n);
else
  I = eye(d^n);
end
end

