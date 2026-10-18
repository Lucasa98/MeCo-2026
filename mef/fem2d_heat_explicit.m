function [PHI_vec,Q_vec] = fem2d_heat_explicit(K,C,F,xnode,icone,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando
% esquema temporal explícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción).
% * C: matriz de masa del sistema (sin escalar por la constante de reacción).
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * icone: matriz de conectividad. Indica los 3 ó 4 nodos que integran el 
%   elemento, recorridos en cualquier orden pero en sentido antihorario. 
%   En caso de elementos triangulares, la cuarta columna siempre es -1.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal crítico para método explícito.

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor 
%   escalar asociado a cada nodo de la malla, y su posición dentro del vector 
%   depende de cómo se especificó cada nodo en xnode. Cada columna representa 
%   una iteración del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas 
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------

    % Nodos Dirichlet: fem2d_heat_dirichlet reemplazó su fila de K por la identidad
    isdir = full(sum(K~=0,2) == 1 & diag(K) == 1);

    % rho*cp*C*(phi^{n+1}-phi^n)/dt + K*phi^n = F
    % En nodos Dirichlet la fila de C pasa a ser a*identidad, así phi = F allí.
    a = dt/(model.rho*model.cp);
    M = C;
    M(isdir,:) = 0;
    M(sub2ind(size(M),find(isdir),find(isdir))) = a;

    phi_n = full(model.PHI_n(:));
    F = full(F);

    PHI_vec = zeros(model.nnodes,model.maxit+1);
    Q_vec = zeros(model.nnodes,2*(model.maxit+1));
    PHI_vec(:,1) = phi_n;
    Q_vec(:,1:2) = fem2d_heat_flux(xnode,icone,model,phi_n);

    nit = 0;
    for n = 1 : model.maxit
        phi = M\(a*(F - K*phi_n) + M*phi_n);

        % Error relativo entre las últimas dos iteraciones
        err = norm(phi-phi_n,2)/norm(phi,2);

        phi_n = phi;
        nit = n;
        PHI_vec(:,n+1) = phi;
        Q_vec(:,2*n+1:2*n+2) = fem2d_heat_flux(xnode,icone,model,phi);

        if err < model.tol
            break; % convergencia a estado estacionario
        end
    end

    PHI_vec = PHI_vec(:,1:nit+1);
    Q_vec = Q_vec(:,1:2*(nit+1));
end