function [K,F] = fdm2d_robin(K,F,xnode,neighb,ROB)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * neighb: matriz de vecindad.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%   - Columna 1: índice del nodo donde se aplica la condición de borde.
%   - Columna 2: valor de coeficiente de calor (h)
%   - Columna 3: valor de temperatura de referencia (phi_inf).
%   - Columna 4: dirección y sentido del flujo:
%     (1) Flujo en dirección eje-y, sentido negativo (S – South – Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East – Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de 
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición 
%   de borde.
% ----------------------------------------------------------------------

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
end