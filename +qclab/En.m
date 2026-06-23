%> @file En.m
%> @brief Implements En with dimension d x d
% ==============================================================================
%> @brief En
%
%> @param n 
%> @param d 
%> @param issparse 
%
% (C) Copyright Daan Camps, Sophia Keip and Roel Van Beeumen 2025.  
% ==============================================================================
function [En] = En(n, d, issparse)
if nargin < 2, d = 2; end
if nargin < 3, issparse = false; end
assert(n < d);
if issparse 
  En = sparse(diag(0:d-1 == n)) ;
else
  En = diag(0:d-1 == n);
end
end

