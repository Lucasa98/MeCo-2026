close all; clear all; more off;

% Problema modelo 1D:
%   rho*cp*dT/dt - k*d2T/dx2 + c*T = G,   0 <= x <= 1,  t >= 0
%   T(0,t) = 0                                  (Dirichlet)
%   dT/dx(1,t) = -T(1,t)  <=>  -k*dT/dn = h*(T - T_inf), h = k = 1, T_inf = 0  (Robin)
%   T(x,0) = 10*x*(1-x)^2                       (condición inicial)
%
% Galerkin con T(x) = sum_j N_j(x)*T_j:
%   rho*cp*C*dT/dt + K*T = F
%   K_ij = int( k*dN_i/dx*dN_j/dx + c*N_i*N_j dx ) + h*N_i*N_j|_rob
%   C_ij = int( N_i*N_j dx )
%   F_i  = int( N_i*G dx ) + h*T_inf*N_i|_rob - q*N_i|_neu + Q*N_i|_pun
%
% Matrices elementales exactas (le = longitud del elemento):
%   Lineal:     K = k/le*[1 -1; -1 1]
%               C = le/6*[2 1; 1 2]
%               F = G*le/2*[1; 1]
%   Cuadrático: K = k/(3*le)*[7 -8 1; -8 16 -8; 1 -8 7]
%               C = le/30*[4 2 -1; 2 16 2; -1 2 4]
%               F = G*le/6*[1; 4; 1]
%
% Esquemas temporales:
%   Explícito:  rho*cp*C*(T^{n+1} - T^n)/dt + K*T^n = F
%   Implícito:  (rho*cp/dt*C + K)*T^{n+1} = F + rho*cp/dt*C*T^n

% Inciso: [1] 2 elementos lineales, estacionario
%         [2] 1 elemento cuadrático, estacionario
%         [3] 2 elementos lineales, implícito con dt = 0.1 (1 iteración)
inciso = 1;

xnode = [
  0.0000000000000000;
  0.5000000000000000;
  1.0000000000000000;
];

% Conectividad: 2 columnas => elemento lineal, 3 columnas => cuadrático
% (nodos del elemento cuadrático ordenados [izquierdo, central, derecho])
if inciso == 2
    icone = [
           1,      2,      3;
    ];
else
    icone = [
           1,      2;
           2,      3;
    ];
end

% nodo, valor de T
DIR = [
       1, 0.0000000000000000;
];

% nodo, q = -k*dT/dn (flujo saliente, n normal exterior)
NEU = [
];

% nodo, h, T_inf  (-k*dT/dn = h*(T - T_inf), n normal exterior)
ROB = [
       3, 1.0000000000000000, 0.0000000000000000;
];

% nodo, fuente puntual Q
PUN = [
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);
model.nelem = size(icone,1);

model.k = 1.0000000000000000;
model.c = 0.0000000000000000;

model.G = [
  100.0000000000000000;
  100.0000000000000000;
];
model.G = model.G(1:model.nelem);

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
if inciso == 3
    model.ts = 1;
else
    model.ts = 2;
end

% Parámetros para esquemas temporales
model.rho = 1.0000000000000000;
model.cp = 1.0000000000000000;
model.dt = 0.1000000000000000;
model.maxit =            1;
model.tol = 1.000000e-05;

% Condición inicial
model.PHI_n = 10*xnode.*(1-xnode).^2;

% Solución exacta para comparar ([] para no graficarla)
% Estacionario: -T'' = 100, T(0) = 0, T'(1) = -T(1)  =>  T = -50x^2 + 75x
if model.ts == 0 || model.ts == 1
    exact = [];
else
    exact = @(x) -50*x.^2 + 75*x;
end

disp('Iniciando el método numérico...');

% Funciones de forma en coordenada natural s en [-1,1]
nen = size(icone,2);
if nen == 2
    shape = @(s) [(1-s)/2, (1+s)/2];
    dshape = @(s) [-1/2, 1/2];
else
    shape = @(s) [s*(s-1)/2, 1-s^2, s*(s+1)/2];
    dshape = @(s) [s-1/2, -2*s, s+1/2];
end

% Cuadratura de Gauss de 3 puntos (exacta hasta grado 5)
gp = [-sqrt(3/5), 0, sqrt(3/5)];
w = [5/9, 8/9, 5/9];

K = zeros(model.nnodes);
C = zeros(model.nnodes);
F = zeros(model.nnodes,1);

% Ensamble de las matrices y vectores del sistema
for e = 1 : model.nelem
    elem = icone(e,:);
    nodes = xnode(elem);

    localK = zeros(nen);
    localC = zeros(nen);
    localF = zeros(nen,1);
    for p = 1 : numel(gp)
        N = shape(gp(p));
        dNds = dshape(gp(p));
        J = dNds*nodes; % dx/ds
        dNdx = dNds/J;
        localK = localK + w(p)*(model.k*(dNdx'*dNdx) + model.c*(N'*N))*abs(J);
        localC = localC + w(p)*(N'*N)*abs(J);
        localF = localF + w(p)*N'*model.G(e)*abs(J);
    end

    K(elem,elem) = K(elem,elem) + localK;
    C(elem,elem) = C(elem,elem) + localC;
    F(elem) = F(elem) + localF;

    fprintf('\nElemento %d, nodos %s:\n', e, mat2str(elem));
    disp('localK ='); disp(rats(localK));
    disp('localC ='); disp(rats(localC));
    disp('localF ='); disp(rats(localF));
end

% Ensamble de nodos frontera Neumann: F_P = F_P - q
for i = 1 : size(NEU,1)
    P = NEU(i,1);
    F(P) = F(P) - NEU(i,2);
end

% Ensamble de nodos frontera Robin: K_PP = K_PP + h, F_P = F_P + h*T_inf
for i = 1 : size(ROB,1)
    P = ROB(i,1);
    K(P,P) = K(P,P) + ROB(i,2);
    F(P) = F(P) + ROB(i,2)*ROB(i,3);
end

% Ensamble de nodos con fuentes puntuales: F_P = F_P + Q
for i = 1 : size(PUN,1)
    P = PUN(i,1);
    F(P) = F(P) + PUN(i,2);
end

disp(' ');
disp('Sistema global antes de Dirichlet (K*T = F):');
disp('K ='); disp(rats(K));
disp('C ='); disp(rats(C));
disp('F ='); disp(rats(F));

% Ensamble de nodos frontera Dirichlet: se reemplaza la ecuación por T_P = valor
isdir = false(model.nnodes,1);
if ~isempty(DIR)
    P = DIR(:,1);
    isdir(P) = true;
    K(P,:) = 0;
    K(sub2ind(size(K),P,P)) = 1;
    F(P) = DIR(:,2);
end

disp('Sistema global con Dirichlet:');
disp('K ='); disp(rats(K));
disp('F ='); disp(rats(F));

% Resolución del sistema lineal de ecuaciones
if model.ts == 0 || model.ts == 1
    disp('Iniciando esquemas temporales...');
    dt = model.dt;
    phi_n = model.PHI_n(:);

    if model.ts == 0
        % Explícito: en nodos Dirichlet la fila de M es a*identidad, así T = F allí
        a = dt/(model.rho*model.cp);
        M = C;
        M(isdir,:) = 0;
        M(sub2ind(size(M),find(isdir),find(isdir))) = a;
    else
        % Implícito: se anulan las filas de M en nodos Dirichlet para que T = F allí
        M = model.rho*model.cp/dt*C;
        M(isdir,:) = 0;
        A = M + K;
        disp('Sistema implícito (rho*cp/dt*C + K)*T^{n+1} = F + rho*cp/dt*C*T^n:');
        disp('A ='); disp(rats(A));
        disp('b (primera iteración) ='); disp(rats(F + M*phi_n));
    end

    PHI = zeros(model.nnodes,model.maxit+1);
    PHI(:,1) = phi_n;

    nit = 0;
    for n = 1 : model.maxit
        if model.ts == 0
            phi = M\(a*(F - K*phi_n) + M*phi_n);
        else
            phi = A\(F + M*phi_n);
        end

        % Error relativo entre las últimas dos iteraciones
        err = norm(phi-phi_n,2)/norm(phi,2);

        phi_n = phi;
        nit = n;
        PHI(:,n+1) = phi;

        if err < model.tol
            break; % convergencia a estado estacionario
        end
    end
    PHI = PHI(:,1:nit+1);
else
    disp('Resolución del sistema de ecuaciones...');
    PHI = K\F;
end

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% Temperatura nodal: columna 1 = x, resto = T en cada paso temporal
disp('    x         T');
disp([xnode, PHI]);

% Interpolación con las funciones de forma dentro de cada elemento
s = linspace(-1,1,21);
figure; hold on; grid on;
for n = 1 : size(PHI,2)
    xs = [];
    Ts = [];
    for e = 1 : model.nelem
        elem = icone(e,:);
        for j = 1 : numel(s)
            N = shape(s(j));
            xs(end+1) = N*xnode(elem);
            Ts(end+1) = N*PHI(elem,n);
        end
    end
    plot(xs,Ts,'-');
    plot(xnode,PHI(:,n),'o');
end
if ~isempty(exact)
    xe = linspace(min(xnode),max(xnode),101);
    plot(xe,exact(xe),'k--');
end
xlabel('x'); ylabel('T');
title(sprintf('MEF 1D - inciso %d', inciso));

disp('Finalizado el post-procesamiento.');
