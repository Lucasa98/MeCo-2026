clear;clc;

f = @(x) [zeros(size(x)); ones(size(x)); zeros(size(x))];
inter = [0 1];
y0 = [0 1];

% h=1/3
[x, y] = diferencias_finitas(f, inter, y0, 2);
[x(2) y(2)]
[x(3) y(3)]

% h=1/6
[x, y] = diferencias_finitas(f, inter, y0, 5);
[x(3) y(3)]
[x(5) y(5)]

