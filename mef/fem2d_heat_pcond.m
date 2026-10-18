function [F] = fem2d_heat_pcond(F,xnode,icone,PUN)
% Descripción: módulo para calcular y ensamblar las contribuciones producidas
% por fuentes puntuales aplicadas en posiciones (x,y) dentro del dominio. Estos
% pares (x,y) también pueden ubicarse en los bordes del elemento, incluso en 
% alguno de sus vértices.

% Entrada:
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * icone: matriz de conectividad. Indica los 3 ó 4 nodos que integran el 
%   elemento, recorridos en cualquier orden pero en sentido antihorario. 
%   En caso de elementos triangulares, la cuarta columna siempre es -1.
% * PUN: matriz con la información sobre las fuentes puntuales aplicadas a 
%   elementos del dominio. 
%   - Columna 1: número de elemento al que se aplica la fuente.
%   - Columna 2: valor de la fuente aplicada (G). 
%   - Columnas 3-4: posición absoluta (x,y) donde se aplica la fuente 
%     (cualquier parte del dominio).

% Salida:
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición
%   de borde.
% ----------------------------------------------------------------------

    for i = 1 : size(PUN,1)
        e = PUN(i,1);
        if icone(e,4) == -1
            elem = icone(e,1:3);
        else
            elem = icone(e,:);
        end
        % Se reparte la fuente en los nodos según las funciones de forma en (x,y)
        N = fem2d_heat_blerp(xnode(elem,:),PUN(i,3),PUN(i,4));
        F(elem) = F(elem) + PUN(i,2)*N;
    end
end
