clc; close all; clear all; more off;

% Problema 2D no estacionario resuelto con fdm2d:
%   rho cp dT/dt - k (d2T/dx2 + d2T/dy2) = G,   0 <= x,y <= 1,  t >= 0
%   k = 10, G = 500, sin reaccion (c = 0), rho cp = 1
%   Borde superior  (y = 1): T = 150                            (Dirichlet)
%   Borde inferior  (y = 0): -k dT/dn = h (T - T_inf),
%                            h = 50, T_inf = 20                 (Robin)
%   Borde izquierdo (x = 0): T = 100                            (Dirichlet)
%   Borde derecho   (x = 1): dT/dx = 0                          (Neumann, aislado)
%   T(x,y,0) = 125
% Malla uniforme dx = dy = 1/6 (7x7 nodos), esquema explicito, tolerancia 1e-5.
% Se informa: tiempo hasta el estado estacionario, temperatura en puntos
% especificos y flujo de calor en el centro de la placa.
%
% Esquinas: donde Dirichlet se toca con Robin o Neumann manda Dirichlet. En (0,1)
% se tocan T = 100 y T = 150 y se usa el promedio. Con el stencil de 5 puntos esa
% esquina no es vecina de ningun nodo libre, asi que su valor no cambia el resto
% de la solucion.

L  = 1;
dx = 1/6;
dy = dx;
nx = round(L/dx) + 1;       % nodos en x
ny = round(L/dy) + 1;       % nodos en y
x  = linspace(0, L, nx);    % valores de x
y  = linspace(0, L, ny);    % valores de y

% Numeracion por columnas (como en main.m): nodo (i-1)*ny + j en (x_i, y_j)
[X, Y] = meshgrid(x, y);
xnode  = [X(:), Y(:)];
id     = reshape(1:nx*ny, ny, nx);   % id(j,i) = nodo en (x_i, y_j)

% Cada elemento: [abajo-izq, abajo-der, arriba-der, arriba-izq]
% 7 - 14 - 21 - ...
% |    |    |   ...
% 2 -  9 - 16 - ...
% |    |    |   ...
% 1 -  8 - 15 - ...
icone = [id(1:end-1, 1:end-1)(:), id(1:end-1, 2:end)(:), ...
         id(2:end,   2:end)(:),   id(2:end,   1:end-1)(:)];

abajo  = id(1, :)';         % nodos en y = 0
arriba = id(end, :)';       % nodos en y = 1
izq    = id(:, 1);          % nodos en x = 0
der    = id(:, end);        % nodos en x = 1

T_izq = 100;
T_sup = 150;
T_esq = (T_izq + T_sup) / 2;   % esquina (0,1)

% DIR [nodo, valor]
DIR = [izq(1:end-1),  T_izq * ones(ny-1, 1);   % incluye la esquina (0,0)
       arriba(2:end), T_sup * ones(nx-1, 1);   % incluye la esquina (1,1)
       id(end, 1),    T_esq];

% NEU [nodo, q, direccion (S, E, N, W)]
% Borde derecho aislado: q = 0 en direccion E (2). Es lo mismo que no imponer
% nada, se deja explicito. Sin la esquina (1,1), que es Dirichlet.
NEU = [der(1:end-1), zeros(ny-1, 1), 2*ones(ny-1, 1)];

% ROB [nodo, h, T_inf, direccion (S, E, N, W)]
% Borde inferior: h = 50, T_inf = 20, flujo en direccion S (1). Sin la esquina
% (0,0), que es Dirichlet.
ROB = [abajo(2:end), 50*ones(nx-1, 1), 20*ones(nx-1, 1), ones(nx-1, 1)];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);

model.k = 10*ones(model.nnodes,1);      % k = 10
model.c = zeros(model.nnodes,1);        % c = 0 (sin reaccion)
model.G = 500*ones(model.nnodes,1);     % G = 500

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 0;

% Parámetros para esquemas temporales
% (en explicito dt lo fija fdm2d_explicit_delta_t, no se usa model.dt)
model.rho = 1;
model.cp = 1;
model.maxit = 10000;
model.tol = 1e-5;

% Condición inicial
% T(x,y,0) = 125
model.PHI_n = 125 * ones(model.nnodes,1);
% fdm2d_explicit no fuerza los nodos Dirichlet: su fila de K es la identidad, asi
% que T^{n+1} = T^n + dt/(rho cp) (T_D - T^n) y solo se relajan hacia el valor de
% contorno con constante de tiempo rho cp = 1, mucho mas lenta que la difusion
% (k / (rho cp) = 10). Se imponen desde t = 0 para que el contorno valga desde
% el inicio y no distorsione el tiempo hasta el estacionario.
model.PHI_n(DIR(:,1)) = DIR(:,2);

disp('Iniciando el método numérico...');

% Llamada principal al Método de Diferencias Finitas
[PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

nit = size(PHI,2) - 1;
dt  = fdm2d_explicit_delta_t(xnode, model);   % el mismo dt que usa fdm2d_solve
t   = (0:nit) * dt;

T  = PHI(:, end);           % estado estacionario
Tg = reshape(T, ny, nx);    % Tg(j,i) = T(x_i, y_j)

printf('dt = %g, pasos hasta tolerancia = %d\n', dt, nit);
printf('Tiempo hasta el estado estacionario: t = %.6f\n', t(end));

% Tabla de temperaturas, con y = 1 arriba como se ve la placa
disp('Temperatura en los nodos en el estado estacionario:');
printf('   y \\ x '); printf('%10.4f', x); printf('\n');
for j = ny:-1:1
    printf('%9.4f ', y(j)); printf('%10.4f', Tg(j,:)); printf('\n');
end

% Temperatura en puntos especificos (bilineal entre nodos, exacta sobre un nodo)
puntos = [
    1/2, 1/2;   % centro
    1/3, 1/3;
    2/3, 2/3;
    1/2, 0;     % medio del borde inferior (Robin)
    1,   1/2;   % medio del borde derecho (aislado)
];
Tp = interp2(X, Y, Tg, puntos(:,1), puntos(:,2));
for p = 1:rows(puntos)
    printf('T(%.4f, %.4f) = %.6f\n', puntos(p,1), puntos(p,2), Tp(p));
end

% Flujo de calor en el centro (fdm2d_flux: diferencias centradas, q = -k grad T)
centro = find(abs(xnode(:,1) - 0.5) < 1e-12 & abs(xnode(:,2) - 0.5) < 1e-12);
qx = Q(centro, end-1);
qy = Q(centro, end);
printf('Flujo de calor en el centro (1/2, 1/2):\n');
printf('  qx = %.6f, qy = %.6f, |q| = %.6f\n', qx, qy, norm([qx, qy]));

figure('Name', 'T(1/2,1/2,t)');
plot(t, PHI(centro, :), 'b.-');
xlabel('t'); ylabel('T(1/2,1/2,t)'); grid on;

% Solo el estado final: con todas las columnas fdm2d_graph_mesh anima cada paso
fdm2d_graph_mesh(T, Q(:, end-1:end), xnode, icone, 0, 0);   % temperatura
fdm2d_graph_mesh(T, Q(:, end-1:end), xnode, icone, 0, 1);   % flujo de calor

disp('Finalizado el post-procesamiento.');

% Error relativo del final del transitorio respecto al estacionario directo
disp('---------------------------------------------------------------');
disp('Calculando el error relativo respecto a la solucion estacionaria...');

model_est    = model;
model_est.ts = 2;   % estacionario
PHI_est = fdm2d(xnode, icone, DIR, NEU, ROB, model_est);

err_2   = norm(T - PHI_est) / norm(PHI_est);
err_inf = norm(T - PHI_est, Inf) / norm(PHI_est, Inf);

printf('Error relativo (norma 2)   = %e\n', err_2);
printf('Error relativo (norma inf) = %e\n', err_inf);
