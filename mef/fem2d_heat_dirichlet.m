function [K,F] = fem2d_heat_dirichlet(K,F,DIR)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos
% pertenecientes a fronteras de tipo Dirichlet.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * DIR: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: número de nodo.
%   - Columna 2: valor en ese nodo (escalar).

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición
%   de borde.
% ----------------------------------------------------------------------

    if isempty(DIR)
        return; % sin fronteras Dirichlet
    end

    P = DIR(:,1);

    % Se reemplaza la ecuación de cada nodo por phi_P = valor
    K(P,:) = 0;
    K(sub2ind(size(K),P,P)) = 1;
    F(P) = DIR(:,2);
end
