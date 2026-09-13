% difFinitasTest.m
% Template de casos de prueba para difFinitas.m
%
% Ecuacion:   rho*c_p * dT/dt = k*T'' - c*T + G(x)     en  L1 < x < L2
%   (si rho_cp = 0 es el estacionario:  k*T'' - c*T + G(x) = 0)
%
% Condiciones de borde en x = L1 (cb1) y x = L2 (cb2), una struct por borde:
%   Dirichlet:  struct("tipo", "dirichlet", "T", T_borde)
%   Neumann:    struct("tipo", "neumann",   "q", q)
%                 q = k*dT/dn con n normal exterior  (q > 0: calor ENTRANDO al dominio)
%                 en L1:  q = -k*T'(L1)      en L2:  q = k*T'(L2)
%   Robin:      struct("tipo", "robin", "h", h, "T_ext", T_ext)
%                 q = h*(T_ext - T_borde)
%
% Campos de cada caso:
%   nombre   : etiqueta
%   L1, L2   : intervalo
%   cb1, cb2 : condicion de borde en L1 y L2
%   rho_cp   : rho*c_p constante. 0 => estacionario.
%              Distinto de 0 => esquema temporal: se avanza hasta el estado
%              estacionario (corte por ET.tol) y se compara contra T.
%   k, c     : constantes
%   G        : fuente, funcion de x (vectorizada)
%   T        : solucion analitica (estacionaria), funcion de x (vectorizada)
%   dt, T0   : (opcionales, solo transitorio) paso temporal y condicion
%              inicial como funcion de x. Si no estan se usa ET.dt y T0 = 0.

clear; clc; close all;

% ================= configuracion =================
Ns       = [10 20 40 80 160];   % nodos interiores de cada refinamiento
graficar = true;            % grafico analitica vs numerica (con el N mas fino)

% esquema temporal (solo para casos con rho_cp != 0)
ET.tipo  = 2;      % 1=Forward Euler, 2=Backward Euler, 3=Crank-Nicolson
ET.maxIt = 1e5;
ET.tol   = 1e-10;  % corte por estado estacionario: ||T^{n+1} - T^n||_inf < tol
ET.dt    = 0;      % 0 => difFinitas usa dt = dtcr/2 (Forward Euler siempre usa 0)

% ================= casos =================
casos = {};

% Item a)
casos{end+1} = struct( ...
    "nombre", "Item a", ...
    "L1", 0, "L2", 1, ...
    "cb1", struct("tipo", "dirichlet", "T", 10), ...
    "cb2", struct("tipo", "dirichlet", "T", 50), ...
    "rho_cp", 0, "k", 2, "c", 0, ...
    "G", @(x) ones(size(x))*100, ...
    "T", @(x) -25*x.^2 + 65*x + 10);

% Item b)
casos{end+1} = struct( ...
    "nombre", "Item b", ...
    "L1", 0, "L2", 2, ...
    "cb1", struct("tipo", "dirichlet", "T", 100), ...
    "cb2", struct("tipo", "neumann", "q", 0), ...
    "rho_cp", 0, "k", 1, "c", 1, ...
    "G", @(x) zeros(size(x)), ...
    "T", @(x) (100.*exp(-x).*(exp(2*x) + exp(4))) / (1 + exp(4)));

% Item c)
casos{end+1} = struct( ...
    "nombre", "Item c", ...
    "L1", 1, "L2", 5, ...
    "cb1", struct("tipo", "neumann", "q", 2), ...
    "cb2", struct("tipo", "dirichlet", "T", 0), ...
    "rho_cp", 0, "k", 1, "c", 0, ...
    "G", @(x) 100*(x-3).^2, ...
    "T", @(x) (-25*x.^4 + 300*x.^3 - 1350*x.^2 + 1906*x + 2345) / 3);

% Item d)
casos{end+1} = struct( ...
    "nombre", "Item d", ...
    "L1", 0, "L2", 1, ...
    "cb1", struct("tipo", "dirichlet", "T", 10), ...
    "cb2", struct("tipo", "robin", "h", 0.2, "T_ext", 50), ...
    "rho_cp", 0, "k", 1, "c", 1, ...
    "G", @(x) ones(size(x))*50, ...
    "T", @(x) -36.6897.*exp(-x) - 3.3103.*exp(x) + 50);

% Item e)
casos{end+1} = struct( ...
    "nombre", "Item e", ...
    "L1", 5, "L2", 10, ...
    "cb1", struct("tipo", "robin", "h", 2, "T_ext", 100), ...
    "cb2", struct("tipo", "dirichlet", "T", 50), ...
    "rho_cp", 1, "k", 2, "c", 0, ...
    "G", @(x) x.^3, ...
    "T", @(x) (-x.^5/40) + (1225*x/3) - (4600/3), ...
    "dt", 1e-3, ...
    "T0", @(x) zeros(size(x)));

% Item f)
casos{end+1} = struct( ...
    "nombre", "Item f", ...
    "L1", 0, "L2", 1, ...
    "cb1", struct("tipo", "dirichlet", "T", 0), ...
    "cb2", struct("tipo", "robin", "h", 2, "T_ext", 10), ...
    "rho_cp", 2, "k", 2, "c", 2, ...
    "G", @(x) ones(size(x))*75, ...
    "T", @(x) -1.25.*exp(-x-1).*(exp(x)-1).*(11.*exp(x) + 11 - 30*e), ...
    "dt", 1e-3, ...
    "T0", @(x) zeros(size(x)));

% Item g)
casos{end+1} = struct( ...
    "nombre", "Item g", ...
    "L1", 0, "L2", 1, ...
    "cb1", struct("tipo", "dirichlet", "T", 50), ...
    "cb2", struct("tipo", "neumann", "q", 5), ...
    "rho_cp", 1, "k", 2, "c", -2, ...
    "G", @(x) zeros(size(x)), ...
    "T", @(x) 73.2433*sin(x) + 50*cos(x), ...
    "dt", 1e-3, ...
    "T0", @(x) zeros(size(x)));

% ================= no hace falta tocar de aca para abajo =================

function fila = armarCB(cb)
    % struct de condicion de borde -> fila de la matriz cb de difFinitas
    switch lower(cb.tipo)
        case "dirichlet"
            fila = [1, cb.T, -1];
        case "neumann"
            fila = [2, cb.q, -1];
        case "robin"
            fila = [3, cb.h, cb.T_ext];
        otherwise
            error("difFinitasTest: tipo de condicion de borde desconocido: %s", cb.tipo);
    endswitch
endfunction

function [T, xnode] = correrCaso(caso, N, ET)
    xnode = linspace(caso.L1, caso.L2, N + 2);
    model = struct("k", caso.k, "c", caso.c, "rho", caso.rho_cp, "c_p", 1, "G", caso.G);
    cb    = [armarCB(caso.cb1); armarCB(caso.cb2)];

    if (caso.rho_cp == 0)
        T = difFinitas(xnode, model, cb, 0);
    else
        dt = ET.dt;
        if (isfield(caso, "dt")), dt = caso.dt; endif
        if (ET.tipo == 1), dt = 0; endif   % Forward Euler: siempre dt = dtcr/2 (estable)
        T0 = [];
        if (isfield(caso, "T0")), T0 = caso.T0(xnode(:)); endif
        T = difFinitas(xnode, model, cb, [ET.tipo, ET.maxIt, ET.tol, dt], T0);
    endif
endfunction

for i = 1:numel(casos)
    caso = casos{i};
    printf("\n=== Caso %d: %s ===\n", i, caso.nombre);
    if (caso.rho_cp == 0)
        printf("estacionario\n");
    else
        printf("transitorio: rho*c_p = %g, esquema %d\n", caso.rho_cp, ET.tipo);
    endif
    printf("   N        h       error max     orden\n");
    errPrev = NaN; hPrev = NaN;
    for N = Ns
        [T, xnode] = correrCaso(caso, N, ET);
        h   = (caso.L2 - caso.L1) / (N + 1);
        err = max(abs(T(:) - caso.T(xnode(:))));
        if (err < 1e-12)
            printf("%4d  %9.5f   %.4e   (DF exacta)\n", N, h, err);
        elseif (isnan(errPrev))
            printf("%4d  %9.5f   %.4e\n", N, h, err);
        else
            printf("%4d  %9.5f   %.4e   %6.2f\n", N, h, err, log(errPrev/err)/log(hPrev/h));
        endif
        errPrev = err; hPrev = h;
    endfor

    if (graficar)
        % volver a correr el primero, que es el que quiero graficar
        [T, xnode] = correrCaso(caso, Ns(1), ET);
        figure(i);
        xf = linspace(caso.L1, caso.L2, 400);
        plot(xf, caso.T(xf), "b-", xnode, T, "ro");
        legend("analitica", "difFinitas");
        title(sprintf("Caso %d: %s", i, caso.nombre));
        xlabel("x"); ylabel("T"); grid on;
    endif
endfor
