function [Qn] = fem2d_heat_flux(xnode,icone,model,PHI)
% Descripción: módulo calcular el flujo de calor en todo el dominio. Se aplica
% la Ley de Fourier y se evalúa como fluye el calor en todos los puntos (nodos)
% del dominio.

% Entrada:
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de 
%   cada nodo de la malla.
% * icone: matriz de conectividad (la cuarta columna es -1 en triángulos).
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * PHI: vector solución. Cada elemento del vector representa un valor escalar
%   asociado a cada nodo de la malla, y su posición dentro del vector depende 
%   de cómo se especificó cada nodo en xnode.

% Salida:
% * Qn: matriz [nnodes×2] de flujo de calor nodal (Qx,Qy).
% ----------------------------------------------------------------------

    phi = full(PHI(:));
    grad = zeros(model.nnodes,2);
    cnt = zeros(model.nnodes,1);

    for e = 1 : model.nelem
        if icone(e,4) == -1
            elem = icone(e,1:3);
            % Triángulo: gradiente constante en el elemento
            dN = [-1 1 0; -1 0 1];
            B = (dN*xnode(elem,:))\dN;
            g = (B*phi(elem))';
            grad(elem,:) = grad(elem,:) + repmat(g,3,1);
        else
            elem = icone(e,:);
            % Cuadrángulo: gradiente evaluado en cada vértice (coord. naturales)
            st = [-1 -1; 1 -1; 1 1; -1 1];
            for a = 1 : 4
                s = st(a,1); t = st(a,2);
                dN = 0.25*[-(1-t) (1-t) (1+t) -(1+t);
                           -(1-s) -(1+s) (1+s) (1-s)];
                B = (dN*xnode(elem,:))\dN;
                grad(elem(a),:) = grad(elem(a),:) + (B*phi(elem))';
            end
        end
        cnt(elem) = cnt(elem) + 1;
    end

    % Promedio nodal del gradiente y Ley de Fourier: q = -k grad(phi)
    grad = grad./cnt;
    Qn = [-model.kx*grad(:,1), -model.ky*grad(:,2)];
end

