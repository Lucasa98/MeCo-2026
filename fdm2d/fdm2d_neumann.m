function [F] = fdm2d_neumann(F, xnode, neighb, NEU)
    % sumar contribuciones de nodos en fronteras Neumann
    % F: vector de flujo termico
    % xnode: matriz de nodos con pares de coordenadas (x,y)
    % neighb: matriz de vecindad
    % NEU: matriz con info de la frontera Neumann
    %   - col 1: indice del nodo donde aplica
    %   - col 2: flujo q
    %   - col 3: direccion y sentido del flujo (negativo: S o W, positivo, E o N)

    if (isempty(NEU))
        return;   % sin fronteras Neumann
    endif

    % el nodo opuesto de cada vecino, para S, N, etc
    OPUESTO = [3 4 1 2];   % S<->N, E<->W
    EJE     = [2 1 2 1];   % eje de cada vecino

    P   = NEU(:, 1);
    q   = NEU(:, 2);
    dir = NEU(:, 3);

    dir_op = OPUESTO(dir)(:);
    eje    = EJE(dir)(:);
    % ni idea, pero la dist al ficticio es la misma que al opuesto real
    op     = neighb(sub2ind(size(neighb), P, dir_op));
    d      = abs(xnode(sub2ind(size(xnode), P, eje)) - xnode(sub2ind(size(xnode), op, eje)));

    % tampoco entiendo, pero un nodo esquina aparece dos veces y sparse acumula
    F = F + sparse(P, 1, -2 * q ./ d, rows(xnode), 1);
endfunction
