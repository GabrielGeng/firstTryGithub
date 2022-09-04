function sa = square_array(J, dx, dy )
%SQUARE_ARRAY(J, dx, dy) Generate a 2D array with four elements
%   Generate a square array with:
%     - sensor spacing in the x direction of dx
%     - sensor spacing in the y direction of dy
%     - center of the array in the origin

% Calculate the sensor positions in a matrix with all x-axis locations in
% the first column and all y-axis locations in the second column.

%p = [px, py];

% Generate the square array with sensor positions p
p = zeros(J,2);
px = (-J/8:J/8).'*dx;    %size J/4 +1
py = (-J/8:J/8).'*dy;
p(:,1) = [px;px];
p(:,end) = [py(2);py(2);py(1);py(1)];
sa = array(p, 'Square array');
end

