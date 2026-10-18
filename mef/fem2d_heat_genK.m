function [localK] = fem2d_heat_genK(nodes,kx,ky)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de
% difusión K. Se utilizan funciones de forma en coordenadas naturales y se 
% resuelve la integral de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * kx: conductividad térmica orientada en eje-x.
% * ky: conductividad térmica orientada en eje-y.

% Salida:
% * localK: matriz de difusión del elemento (local).
% ----------------------------------------------------------------------

    k = [kx 0; 0 ky];

    if size(nodes,1) == 3
        % Triángulo: B y J constantes, área del triángulo normalizado = 1/2
        dN = [-1 1 0; -1 0 1];
        J = dN*nodes;
        B = J\dN;
        localK = B'*k*B*abs(det(J))/2;
    else
        % Cuadrángulo: cuadratura de Gauss 2x2 con pesos w = 1
        g = 1/sqrt(3);
        gp = [-g -g; g -g; -g g; g g];
        localK = zeros(4);
        for p = 1 : 4
            s = gp(p,1); t = gp(p,2);
            dN = 0.25*[-(1-t) (1-t) (1+t) -(1+t);
                       -(1-s) -(1+s) (1+s) (1-s)];
            J = dN*nodes;
            B = J\dN;
            localK = localK + B'*k*B*abs(det(J));
        end
    end
end
