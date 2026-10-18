function [F] = fdm2d_neumann(F,xnode,neighb,NEU)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Neumann.

% Entrada:
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada 
%   nodo de la malla.
% * neighb: matriz de vecindad.
% * NEU: matriz con la información sobre la frontera de tipo Neumann. 
%   - Columna 1: índice del nodo donde se aplica la condición de borde.
%   - Columna 2: valor de flujo térmico (q) asociado al lado del elemento.
%   - Columna 3: dirección y sentido del flujo:
%     (1) Flujo en dirección eje-y, sentido negativo (S – South - Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East - Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North - Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición de borde.
% ----------------------------------------------------------------------

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
end