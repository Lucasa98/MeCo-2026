function [PHI, Q] = fdm2d_explicit(K, F, xnode, neighb, model, dt)
    % Resolucion del sistema con esquema temporal explicito (Forward Euler), o sea
    % la ec. (168) con theta = 0. Despejando se llega a la ec. (170):
    %
    %   phi^{n+1} = A*F + (I - A*K)*phi^n      con A = dt/(rho*cp)
    %
    % Es condicionalmente estable: dt debe respetar el paso critico que calcula
    % fdm2d_explicit_delta_t. La primer columna de la salida es la condicion
    % inicial, y su par de columnas en Q es cero (aun no se evaluo ningun flujo).
    %
    % K      [nnodesxnnodes sparse double]: matriz del sistema (difusion + reaccion)
    % F      [nnodesx1      sparse double]: vector de flujo termico
    % xnode  [nnodesx2      double]: coordenadas (x,y), para evaluar el flujo
    % neighb [nnodesx4      int]: vecindad [S E N W], para evaluar el flujo
    % model  [1x1           struct]: se usan nnodes, k, rho, cp, maxit, tol y PHI_n
    % dt     [1x1           double]: paso temporal (critico, del metodo explicito)
    % PHI    [nnodesx(nit+1)     double]: una columna por iteracion
    % Q      [nnodesx(2*(nit+1)) double]: un par (Qx,Qy) por columna de PHI

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
endfunction
