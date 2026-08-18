% Prueba de diferencias_finitas.m
%
%   y'' = -(2/x) y' + (2/x^2) y + sin(ln x)/x^2 ,  1 <= x <= 2
%   y(1) = 1,  y(2) = 2
%
% Solucion exacta:
%   y(x) = c1*x + c2/x^2 - (3/10)sin(ln x) - (1/10)cos(ln x)

clear; clc;

f     = @(x) [-2./x ; 2./x.^2 ; sin(log(x))./x.^2];
inter = [1 2];
y0    = [1 2];

c1 = 1.1392070132;
c2 = -0.03920701320;
exacta = @(x) c1*x + c2*x.^(-2) - 0.3*sin(log(x)) - 0.1*cos(log(x));

printf("  N       h        error maximo    razon\n");
prev = NaN;
for N = [9 19 39 79 159]
  [x, y] = diferencias_finitas(f, inter, y0, N);
  err = max(abs(y - exacta(x)));
  if isnan(prev)
    printf("%3d  %8.5f   %.6e\n", N, (inter(2)-inter(1))/(N+1), err);
  else
    printf("%3d  %8.5f   %.6e   %6.3f\n", N, (inter(2)-inter(1))/(N+1), err, prev/err);
  end
  prev = err;
end

% grafico
[x, y] = diferencias_finitas(f, inter, y0, 19);
xf = linspace(inter(1), inter(2), 400);
plot(xf, exacta(xf), 'b-', x, y, 'ro');
legend('exacta', 'diferencias finitas', 'location', 'northwest');
xlabel('x'); ylabel('y'); grid on;
