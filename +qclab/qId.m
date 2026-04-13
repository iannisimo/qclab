%> @file qId.m
%> @brief Implements quantum Identity matrix
% ==============================================================================
%> @brief Identity matrix on n qubits
%
%> @param n: number of qubits
%> @param issparse: use speye instead of eye (default: false)
%> @param d: number of energy levels (default: 2)
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [I] = qId(n, issparse, d)
if nargin <= 1, issparse = false; end
if nargin <= 2, d = 2; end
if issparse 
  I = speye(d^n) ;
else
  I = eye(d^n);
end
end

