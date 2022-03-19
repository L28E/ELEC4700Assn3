function oob = get_oob (Px, Py, left, right, top, bottom)
  oob = (Px >= left & Px <= right & Py <= top) | (Px >= left & Px <= right & Py >= bottom);
end

