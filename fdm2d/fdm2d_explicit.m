function [PHI, Q] = fdm2d_explicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal explícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo 
% de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal crítico para método explícito.

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor escalar 
%   asociado a cada nodo de la malla, y su posición dentro del vector depende de
%   cómo se especificó cada nodo en xnode. Cada columna representa una iteración
%   del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas 
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------

    nnodes = model.nnodes;
    A = dt / (model.rho * model.cp);

    phi_n = full(model.PHI_n(:));

    PHI = zeros(nnodes, model.maxit + 1);
    Q   = zeros(nnodes, 2 * (model.maxit + 1));
    PHI(:, 1) = phi_n;   % Q(:,1:2) queda en cero: la condicion inicial no aporta flujo

    nit = 0;
    for n = 1:model.maxit
        % equivale a A*F + (I - A*K)*phi^n, sin armar la identidad densa
        phi = A * F + phi_n - A * (K * phi_n);

        % error relativo entre las ultimas dos iteraciones
        err = norm(phi - phi_n, 2) / norm(phi, 2);

        phi_n = phi;
        nit = n;
        PHI(:, n + 1) = phi;
        Q(:, 2*n+1 : 2*n+2) = fdm2d_flux(phi, neighb, xnode, model.k);

        if (err < model.tol)
            break;   % terminado por tolerancia de error
        endif
    endfor

    PHI = PHI(:, 1:nit + 1);
    Q   = Q(:, 1:2 * (nit + 1));
end