function [K,F] = fem2d_heat_robin(K,F,ROB,xnode)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos
% pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%   - Columnas 1-2: dos nodos contiguos formando un lado de un elemento.
%   - Columna 3: valor de coeficiente de calor (h).
%   - Columna 4: valor de temperatura de referencia (phi_inf).
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo de la malla.

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la 
%   condición de borde.
% ----------------------------------------------------------------------

    % -k dphi/dn = h(phi - phi_inf)  =>  K = K + int(h N N dGamma)
    %                                    F = F + int(N h phi_inf dGamma)
    for i = 1 : size(ROB,1)
        n = ROB(i,1:2);
        h = ROB(i,3);
        phi_inf = ROB(i,4);
        d = norm(xnode(n(1),:) - xnode(n(2),:));
        K(n,n) = K(n,n) + h*d/6*[2 1; 1 2];
        F(n) = F(n) + h*phi_inf*d/2;
    end
end
