%> @file isSparse.m
%> @brief Implements isSparse
% ==============================================================================
%> @brief Checks if simulation should be sparse or not dependend on the
%>        number of qubits `nbQubits`
%
%> @param nbQubits number of qubits
%
% (C) Copyright Daan Camps, Sophia Keip and Roel Van Beeumen 2025.  
% ==============================================================================
function [isSparse] = isSparse(nbQubits, d)
  if nargin <= 1, d = 2; end
  if d == 2
    if nbQubits < 9
      isSparse = false ;
    else
      isSparse = true ;
    end
  else
    isSparse = d^nbQubits < 512;
  end
end

