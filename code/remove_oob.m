function [Px, Py] = remove_oob (Px, Py, oob, left, right, top, bottom, M, I)
  Px (oob & M' ~= 0 & I' == 1) = 2 * left - Px (oob & M' ~= 0 & I' == 1);
  Px (oob & M' ~= 0 & I' == 2) = 2 * right - Px (oob & M' ~= 0 & I' == 2);
  Py (oob & M' ~= 0 & I' == 3) = 2 * top - Py (oob & M' ~= 0 & I' == 3);
  Py (oob & M' ~= 0 & I' == 4) = 2 * bottom - Py (oob & M' ~= 0 & I' == 4);
end

