function [K, F] = fdm2d_robin(K, F, xnode, neighb, ROB)
    % calcula contribuciones de fronteras ROBIN
    % K: matriz del sistema
    % F: vector de flujo termico
    % xnode: matriz de nodos con pares (x,y) de coordenadas
    % neighb: matriz de vecindad
    % ROB: matriz con la info de la frontera
    %   - col 1: indice del nodo con frontera Robin
    %   - col 2: coeficiente de calor h
    %   - col 3: temperatura externa phi_inf
    %   - col 4: direccion y sentido del flujo, igual que en fdm2d_neumann

    if (isempty(ROB))
        return;   % el problema no tiene fronteras Robin
    endif

    OPUESTO = [3 4 1 2];   % S<->N, E<->W
    EJE     = [2 1 2 1];   % eje de cada vecino

    P       = ROB(:, 1);
    h       = ROB(:, 2);
    phi_inf = ROB(:, 3);
    dir     = ROB(:, 4);

    dir_op = OPUESTO(dir)(:);
    eje    = EJE(dir)(:);
    op     = neighb(sub2ind(size(neighb), P, dir_op));
    d      = abs(xnode(sub2ind(size(xnode), P, eje)) - xnode(sub2ind(size(xnode), op, eje)));

    nnodes = rows(xnode);
    K = K + sparse(P, P, 2 * h ./ d, nnodes, nnodes);
    F = F + sparse(P, 1, 2 * h .* phi_inf ./ d, nnodes, 1);
endfunction
