function [localC] = fem2d_heat_genC(nodes)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de
% masa C. Se utilizan funciones de forma en coordenadas naturales y se 
% resuelve la integral de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.

% Salida:
% * localC: matriz de masa para elemento (local).
% ----------------------------------------------------------------------

    if size(nodes,1) == 3
        % Triángulo: 3 puntos de Gauss con peso w = 1/6
        gp = [0.5 0; 0 0.5; 0.5 0.5];
        w = 1/6;
        detJ = abs(det([-1 1 0; -1 0 1]*nodes));
        localC = zeros(3);
        for p = 1 : 3
            s = gp(p,1); t = gp(p,2);
            N = [1-s-t s t];
            localC = localC + w*(N'*N)*detJ;
        end
    else
        % Cuadrángulo: cuadratura de Gauss 2x2 con pesos w = 1
        g = 1/sqrt(3);
        gp = [-g -g; g -g; -g g; g g];
        localC = zeros(4);
        for p = 1 : 4
            s = gp(p,1); t = gp(p,2);
            N = 0.25*[(1-s)*(1-t) (1+s)*(1-t) (1+s)*(1+t) (1-s)*(1+t)];
            dN = 0.25*[-(1-t) (1-t) (1+t) -(1+t);
                       -(1-s) -(1+s) (1+s) (1-s)];
            localC = localC + (N'*N)*abs(det(dN*nodes));
        end
    end
end
