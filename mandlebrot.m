%% define step-size
step_size = 0.05;
%% Define the region.
x = 0: step_size: 0.8;
y = x';
%% Create the two-dimensional complex grid
n = length(x);
e = ones(n,1);
z0 = x(e,:) + 1i*y(:,e);
%% You can also do the same thing with meshgrid.
%[X,Y] = meshgrid(x,y);
%z0 = X + 1i*Y;
%% Initialize the iterates and counts arrays.
z = zeros(n,n);
c = zeros(n,n);
%% Mandelbrot iteration.
depth = 32;
for k = 1:depth
    z = z.^3 + z0;
    c(abs(z) < 2) = k;
end
%% Image creaqtion from c
c
image(c)
axis image
%% Colors
colormap(flipud(jet(depth)))