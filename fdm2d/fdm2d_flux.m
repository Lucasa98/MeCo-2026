function [Q] = fdm2d_flux(PHI, neighb, xnode, k)
    % Calculo del flujo de calor en todo el dominio por la Ley de Fourier,
    %
    %   Q = -k * grad(phi)     ->   Qx = -k*dphi/dx,  Qy = -k*dphi/dy
    %
    % Cada derivada se aproxima por diferencias sobre el propio eje: si el nodo
    % tiene sus dos vecinos se usa la diferencia centrada entre ellos, y si es un
    % nodo de frontera se usa la diferencia de dos puntos contra el unico vecino
    % que existe sobre ese eje.
    %
    % PHI    [nnodesx1 sparse double]: solucion sobre la que se evalua el flujo
    % neighb [nnodesx4 int]: vecindad [S E N W], -1 si el vecino no existe
    % xnode  [nnodesx2 double]: coordenadas (x,y) de cada nodo
    % k      [nnodesx1 double]: conductividad termica, distribucion k(x,y)
    % Q      [nnodesx2 double]: par (Qx,Qy) por nodo

    S = 1; E = 2; N = 3; W = 4;

    phi = full(PHI(:));
    k   = k(:);

    dphidx = derivada_eje(phi, xnode(:, 1), neighb, E, W);
    dphidy = derivada_eje(phi, xnode(:, 2), neighb, N, S);

    Q = [-k .* dphidx, -k .* dphidy];
endfunction

function [d] = derivada_eje(phi, x, neighb, pos, neg)
    % Derivada de phi respecto del eje cuyas coordenadas son x, en todos los nodos.
    nnodes = numel(phi);
    vp = neighb(:, pos); vn = neighb(:, neg);
    hayP = (vp != -1); hayN = (vn != -1);

    % Nodo con ambos vecinos: diferencia centrada entre ellos.
    % Nodo de frontera: diferencia entre el nodo y su unico vecino sobre el eje.
    A = zeros(nnodes, 1); B = zeros(nnodes, 1);
    ambos = hayP & hayN;
    A(ambos) = vp(ambos);        B(ambos) = vn(ambos);
    soloP = hayP & !hayN;
    A(soloP) = vp(soloP);        B(soloP) = find(soloP);
    soloN = hayN & !hayP;
    A(soloN) = find(soloN);      B(soloN) = vn(soloN);

    d = (phi(A) - phi(B)) ./ (x(A) - x(B));
endfunction
