function [K, F] = fdm2d_gen_system(K, F, xnode, neighb, k, c, G)
    % Ensamble de los terminos difusivo, reactivo y fuente para todos los nodos de
    % la malla, generando el stencil adecuado segun el nodo sea interior o de
    % frontera. Por cada eje se arman los tres coeficientes de la derivada segunda:
    % 'a' para el vecino positivo (E o N), 'c' para el negativo (W o S) y 'b' para
    % el propio nodo.
    %
    % En un nodo interior salen de la ec. (35), valida para espaciado no uniforme.
    % En un nodo de frontera el vecino que falta es un nodo ficticio ubicado a la
    % misma distancia que su opuesto, con lo que su coeficiente se pliega sobre el
    % vecino real (que duplica su contribucion) y desaparece de K: es el stencil de
    % la ec. (38). La ecuacion que falta la aportan neumann/robin/dirichlet.
    %
    % Los coeficientes se ASIGNAN, no se acumulan: este modulo define la fila de
    % cada nodo, y las demas posiciones de K quedan como estaban.
    %
    % K      [nnodesxnnodes sparse double]: matriz del sistema (difusion + reaccion)
    % F      [nnodesx1      sparse double]: vector de flujo termico
    % xnode  [nnodesx2      double]: coordenadas (x,y) de cada nodo
    % neighb [nnodesx4      int]: vecindad [S E N W], -1 si el vecino no existe
    % k      [nnodesx1      double]: conductividad termica, distribucion k(x,y)
    % c      [nnodesx1      double]: constante de reaccion, distribucion c(x,y)
    % G      [nnodesx1      double]: fuente de calor, distribucion G(x,y)

    S = 1; E = 2; N = 3; W = 4;

    nnodes = size(xnode, 1);
    k = k(:); c = c(:); G = G(:);

    vS = neighb(:, S); vE = neighb(:, E); vN = neighb(:, N); vW = neighb(:, W);
    hayS = (vS != -1); hayE = (vE != -1); hayN = (vN != -1); hayW = (vW != -1);

    % Distancia a cada vecino (0 si no existe)
    ds = zeros(nnodes, 1); de = zeros(nnodes, 1);
    dn = zeros(nnodes, 1); dw = zeros(nnodes, 1);
    ds(hayS) = abs(xnode(vS(hayS), 2) - xnode(hayS, 2));
    de(hayE) = abs(xnode(vE(hayE), 1) - xnode(hayE, 1));
    dn(hayN) = abs(xnode(vN(hayN), 2) - xnode(hayN, 2));
    dw(hayW) = abs(xnode(hayW, 1) - xnode(vW(hayW), 1));

    [ax, bx, cx] = coeficientes(hayE, hayW, de, dw);   % eje x
    [ay, by, cy] = coeficientes(hayN, hayS, dn, ds);   % eje y

    % Todo nodo aporta su reaccion c, los terminos centrales 'bx' y 'by', y su
    % fuente G. Los vecinos solo aportan si existen.
    P = (1:nnodes)';
    K(sub2ind(size(K), P, P)) = c - k .* bx - k .* by;
    F(P) = G;

    K(sub2ind(size(K), P(hayS), vS(hayS))) = -k(hayS) .* cy(hayS);
    K(sub2ind(size(K), P(hayE), vE(hayE))) = -k(hayE) .* ax(hayE);
    K(sub2ind(size(K), P(hayN), vN(hayN))) = -k(hayN) .* ay(hayN);
    K(sub2ind(size(K), P(hayW), vW(hayW))) = -k(hayW) .* cx(hayW);
endfunction

function [a, b, cc] = coeficientes(hayPos, hayNeg, dpos, dneg)
    % Coeficientes de la derivada segunda sobre un eje. Si falta el vecino
    % negativo, el nodo ficticio se pliega sobre el positivo (y viceversa), que es
    % lo que duplica el coeficiente y anula el del ausente.
    n = numel(dpos);
    a = zeros(n, 1); b = zeros(n, 1); cc = zeros(n, 1);

    soloPos = hayPos & !hayNeg;
    a(soloPos)  =  2 ./ (dpos(soloPos) .* dpos(soloPos));
    b(soloPos)  = -2 ./ (dpos(soloPos) .* dpos(soloPos));
    cc(soloPos) =  0;

    soloNeg = hayNeg & !hayPos;
    a(soloNeg)  =  0;
    b(soloNeg)  = -2 ./ (dneg(soloNeg) .* dneg(soloNeg));
    cc(soloNeg) =  2 ./ (dneg(soloNeg) .* dneg(soloNeg));

    amb = !soloPos & !soloNeg;
    a(amb)  =  2 ./ (dpos(amb) .* (dpos(amb) + dneg(amb)));
    b(amb)  = -2 ./ (dpos(amb) .* dneg(amb));
    cc(amb) =  2 ./ (dneg(amb) .* (dpos(amb) + dneg(amb)));
endfunction
