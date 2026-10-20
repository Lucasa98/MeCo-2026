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
%  izq', 0;
];

% ROB [nodo, q, direccion (S, E, N, W)]
NEU = [
  izq', 10;
];

% ROB [nodo, h, T_inf, direccion (S, E, N, W)]
% dT/dx(1,t) = -T(1,t): h = 1, phi_inf = 0, flujo en direccion E (2)
ROB = [
  der', 5, 0;
];

PUN = [
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.kx = 2.0000000000000000;
model.ky = 2.0000000000000000;
model.c = 3.0000000000000000;
model.G = 60*ones(model.nelem,1);

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
% Inciso 1: model.ts = 2;  Inciso 3: model.ts = 1;
model.ts = 2;

% Parámetros para esquemas temporales
model.rho = 2.0000000000000000;
model.cp = 2.0000000000000000;
model.dt = 0.1000000000000000;
model.maxit =            1;
model.tol = 1.000000e-05;

% Condición inicial
model.PHI_n = 10*xnode(:,1).*(1-xnode(:,1)).^2;

disp('Iniciando el método numérico...');

% Llamada principal al Método de Elementos Finitos
[PHI,Q] = fem2d_heat(xnode,icone,DIR,NEU,ROB,PUN,model);

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% Temperatura en los nodos de y = 0: columna 1 = x, resto = T en cada paso
ib = xnode(:,2) == 0;
disp('    x         T');
disp([xnode(ib,1), full(PHI(ib,:))]);

% mode ---> modo de visualización:
%           [0] 2D - Con malla
%           [1] 3D - Con malla
%           [2] 2D - Sin malla
%           [3] 3D - Sin malla
% graph --> tipo de gráfica:
%           [0] Temperatura (escalar)
%           [1] Flujo de Calor (vectorial)
%           [2] Flujo de Calor eje-x (escalar)
%           [3] Flujo de Calor eje-y (escalar)
%           [4] Magnitud de Flujo de Calor (escalar)
mode = 0;
graph = 0;
fem2d_heat_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

disp('Finalizado el post-procesamiento.');
