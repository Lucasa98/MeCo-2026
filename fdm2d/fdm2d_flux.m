function [Q] = fdm2d_flux(PHI,neighb,xnode,k)
% Descripción: módulo calcular el flujo de calor en todo el dominio. Se aplica la
% Ley de Fourier y se evalúa como fluye el calor en todos los puntos (nodos) del dominio.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * PHI: vector solución. Cada elemento del vector representa un valor escalar asociado
% a cada nodo de la malla, y su posición dentro del vector depende de cómo se especificó
% cada nodo en xnode.

% Salida:
% * Q: vector de flujo de calor. Para cada nodo se halla un vector bidimensional de 
%   flujo de calor, representado por un par (Qx,Qy)
% ----------------------------------------------------------------------

    S = 1; E = 2; N = 3; W = 4;

    phi = full(PHI(:));
    k   = k(:);

    dphidx = derivada_eje(phi, xnode(:, 1), neighb, E, W);
    dphidy = derivada_eje(phi, xnode(:, 2), neighb, N, S);

    Q = [-k .* dphidx, -k .* dphidy];
end

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
end
