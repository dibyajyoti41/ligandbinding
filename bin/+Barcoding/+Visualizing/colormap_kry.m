function h = colormap_kry(m)
%colormap_kry    blacK-Red-Yellow color map

if nargin < 1, m = size(get(gcf,'colormap'),1); end
o = m;
m = ceil(m*4/3);
n = fix(3/8*m);

r = [(1:n)'/n; ones(m-n,1)];
g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
b = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];

h = [r g b];

h = h(1:o, :);