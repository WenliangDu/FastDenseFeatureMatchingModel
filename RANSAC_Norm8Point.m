function f = RANSAC_Norm8Point(Location1, Location2)
% Normalize the points
num = size(Location1, 2);
[Location1, Trans1] = vision.internal.normalizePoints(Location1, 2, 'double');
[Location2, Trans2] = vision.internal.normalizePoints(Location2, 2, 'double');

% Compute the constraint matrix
A = double(zeros(num, 9));
for idx = 1: num
    A(idx,:) = [...
        Location1(1,idx)*Location2(1,idx), Location1(2,idx)*Location2(1,idx), Location2(1,idx), ...
        Location1(1,idx)*Location2(2,idx), Location1(2,idx)*Location2(2,idx), Location2(2,idx), ...
        Location1(1,idx),              Location1(2,idx), 1];
end

% Find out the eigen-vector corresponding to the smallest eigen-value.
[~, ~, vA] = svd(A, 0);
f = reshape(vA(:, end), 3, 3)';

% Enforce rank-2 constraint
[u, s, v] = svd(f);
s(end) = 0;
f = u * s * v';

% Transform the fundamental matrix back to its original scale.
f = Trans2' * f * Trans1;

% Normalize the fundamental matrix.
f = f / norm(f);
if f(end) < 0
    f = -f;
end