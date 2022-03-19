function [minimum, index] = get_nearest_bound (num_electrons, Px, Py, oob, left, right, top, bottom)
  diffs = zeros (num_electrons, 4);
  diffs (oob, :) = [abs(Px (oob)' - left), abs(Px (oob)' - right), abs(Py (oob)' - top), abs(Py (oob)' - bottom)];
  [minimum, index] = min (diffs, [], 2);
end

