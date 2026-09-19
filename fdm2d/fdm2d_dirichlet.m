function [K, F] = fdm2d_dirichlet(K, F, DIR)
    % calcular contribuciones de fronteras dirichlet
    % K: matriz del sistema
    % F: vector de flujo termico
    % DIR: matriz con info de la frontera Dirichlet
    %   - col 1: indice del nodo
    %   - col 2: valor de temperatura

    if (isempty(DIR))
        return;   % el problema no tiene fronteras Dirichlet
    endif

    P   = DIR(:, 1);
    val = DIR(:, 2);

    % no acumula, pisa y a la bosta
    K(P, :) = 0;
    K(sub2ind(size(K), P, P)) = 1;
    F(P) = val;
endfunction
