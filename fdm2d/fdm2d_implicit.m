function [PHI, Q] = fdm2d_implicit(K, F, xnode, neighb, model, dt)
    % Resolucion del sistema con esquema temporal implicito (Backward Euler), o sea
    % la ec. (168) con theta = 1. Despejando se llega a la ec. (171):
    %
    %   K_imp = (rho*cp/dt)*I + K        F_imp = F + (rho*cp/dt)*phi^n
    %   K_imp * phi^{n+1} = F_imp
    %
    % es decir un sistema lineal por iteracion. Es incondicionalmente estable, por
    % lo que dt es arbitrario. La primer columna de la salida es la condicion
    % inicial, y su par de columnas en Q es cero (aun no se evaluo ningun flujo).
    %
    % K      [nnodesxnnodes sparse double]: matriz del sistema (difusion + reaccion)
    % F      [nnodesx1      sparse double]: vector de flujo termico
    % xnode  [nnodesx2      double]: coordenadas (x,y), para evaluar el flujo
    % neighb [nnodesx4      int]: vecindad [S E N W], para evaluar el flujo
    % model  [1x1           struct]: se usan nnodes, k, rho, cp, maxit, tol y PHI_n
    % dt     [1x1           double]: paso temporal (arbitrario)
    % PHI    [nnodesx(nit+1)     double]: una columna por iteracion
    % Q      [nnodesx(2*(nit+1)) double]: un par (Qx,Qy) por columna de PHI

    nnodes = model.nnodes;
    a = model.rho * model.cp / dt;

    % speye en vez de eye: si K es sparse el sistema se mantiene sparse y se
    % resuelve como tal, y si K es densa el resultado es denso igual.
    K_imp = K + a * speye(nnodes);

    phi_n = full(model.PHI_n(:));

    PHI = zeros(nnodes, model.maxit + 1);
    Q   = zeros(nnodes, 2 * (model.maxit + 1));
    PHI(:, 1) = phi_n;   % Q(:,1:2) queda en cero: la condicion inicial no aporta flujo

    nit = 0;
    for n = 1:model.maxit
        F_imp = F + a * phi_n;
        phi   = K_imp \ F_imp;

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
