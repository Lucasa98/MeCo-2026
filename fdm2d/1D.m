clc; close all; clear all; more off;

nx = 21;                    % numero de nodos
L  = 1;
x  = linspace(0, L, nx)';   % valores de x
dx = L / (nx - 1);          % delta x
dy = dx;   % arbitrario, no cambia la solucion

% Numeracion por columnas (como en main.m): nodo 2i-1 en (x_i,0), nodo 2i en (x_i,dy)
xnode = zeros(2*nx, 2);
xnode(1:2:end, 1) = x;  % nodos de abajo
xnode(2:2:end, 1) = x;  % nodos de arriba
xnode(2:2:end, 2) = dy; % nodos de arriba

% Cada elemento: [abajo-izq, abajo-der, arriba-der, arriba-izq]
% 2 - 4 - 6 - 8 - ...
% |   |   |   |   ...
% 1 - 3 - 5 - 7 - ...
j = (1:nx-1)';
icone = [2*j-1, 2*j+1, 2*j+2, 2*j];

izq = [1; 2];               % nodos en x = 0
der = [2*nx-1; 2*nx];       % nodos en x = 1

% T(0,t) = 0
DIR = [
%  izq, zeros(2,1);
];

% ROB [nodo, q, direccion (S, E, N, W)]
NEU = [
  izq, 10*ones(2,1), 4*ones(2,1);
];

% ROB [nodo, h, T_inf, direccion (S, E, N, W)]
% dT/dx(1,t) = -T(1,t): h = 1, phi_inf = 0, flujo en direccion E (2)
ROB = [
  der, 5*ones(2,1), zeros(2,1), 2*ones(2,1)
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);

model.k = 2*ones(model.nnodes,1);
model.c = 3*ones(model.nnodes,1);
model.G = 60*ones(model.nnodes,1);

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 2;

% Parámetros para esquemas temporales
model.rho = 2;
model.cp = 2;
lambda = 40;                                            % lambda = k dt / (rho cp dx^2)
model.dt = lambda * model.rho * model.cp * dx^2 / min(model.k);
model.maxit = 100000;
model.tol = 1e-4;

% Condición inicial
% T(x,0) = 10x (1-x)^2
model.PHI_n = 10 * xnode(:,1) .* (1 - xnode(:,1)).^2;

disp('Iniciando el método numérico...');

% Llamada principal al Método de Diferencias Finitas
[PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% Temperatura en los nodos de y = 0: columna 1 = x, resto = T en cada paso
ib = xnode(:,2) == 0;
disp('    x         T');
disp([xnode(ib,1), full(PHI(ib,:))]);

nit = size(PHI,2) - 1;
t   = (0:nit) * model.dt;

% Nodo en x = 1/2 (fila de abajo; la de arriba es identica)
nodo  = find(abs(xnode(:,1) - 0.5) < 1e-12 & xnode(:,2) == 0);
T_med = PHI(nodo, :);

printf('dt = %g, pasos hasta tolerancia = %d (t = %g)\n', model.dt, nit, t(end));
printf('T(1/2) final = %.6f (estacionario exacto: 25)\n', T_med(end));
printf('Diferencia entre filas y=0 e y=dy: %g\n', max(max(abs(PHI(1:2:end,:) - PHI(2:2:end,:)))));

n = find(T_med > 23, 1);
if isempty(n)
    disp('T(1/2,t) no supera 23 antes de alcanzar la tolerancia.');
    n = nit + 1;   % para los graficos
else
    printf('Primer instante con T(1/2,t) > 23: paso %d, t = %g\n', n-1, t(n));
    printf('  T(1/2, %g) = %.6f\n', t(n-1), T_med(n-1));
    printf('  T(1/2, %g) = %.6f\n', t(n), T_med(n));
end

figure('Name', 'T(1/2,t)');
plot(t, T_med, 'b.-', [t(1) t(end)], [23 23], 'r--');
xlabel('t'); ylabel('T(1/2,t)'); grid on;
legend('T(1/2,t)', 'T = 23', 'location', 'southeast');

figure('Name', 'Perfiles');
plot(x, PHI(1:2:end, 1), 'k-', x, PHI(1:2:end, n), 'b-', x, PHI(1:2:end, end), 'r-', ...
     x, 75*x - 50*x.^2, 'r:');
xlabel('x'); ylabel('T'); grid on;
legend('t = 0', sprintf('t = %g', t(n)), sprintf('t = %g', t(end)), 'estacionario exacto', ...
       'location', 'southeast');

disp('Finalizado el post-procesamiento.');

% Error relativo de la solucion estacionaria respecto a la exacta
% T_ex(x) = 25 x (3 - 2x)
disp('---------------------------------------------------------------');
disp('Calculando el error relativo de la solucion estacionaria...');

model_est    = model;
model_est.ts = 2;   % estacionario
PHI_est = fdm2d(xnode, icone, DIR, NEU, ROB, model_est);

% Para comparar con el estado final del transitorio: PHI(1:2:end, end)
T_est = PHI_est(1:2:end);   % fila de abajo.
T_ex  = 25 * x .* (3 - 2*x);

err_2   = norm(T_est - T_ex) / norm(T_ex);
err_inf = norm(T_est - T_ex, Inf) / norm(T_ex, Inf);

printf('Error relativo estacionario (norma 2)   = %e\n', err_2);
printf('Error relativo estacionario (norma inf) = %e\n', err_inf);
