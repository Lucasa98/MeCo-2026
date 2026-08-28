function [T] = difFinitas(xnode, model, cb, et)
  % DIFERENCIAS FINITAS 1D MALLA UNIFORME

  % xnode: vector de coordenadas nodales
  % model: struct con las constantes del modelo (k, c, rho, c_p, G(x))
  % cb: matriz de 2X3
  %   fila 1: condicion de borde del lado izquierdo
  %   fila 2: condicion de borde del lado derecho
  %   columna 1: tipo de condicion (1=Dirichlet, 2=Neumann, 3=Robin)
  %   columna 2: temperatura, flujo o coef de conveccion h (dependiendo de la condicion de la columan 1)
  %   columna 3: -1 para Dirichlet y Neumann, valor de temperatura externa para condicion Robin.
  % et: escalar que indica esquema temporal a utilizar (o si se resuelve en estado estacionario)

  % Casos sin solucion (sin Dirichlet en ningun lado)
  if (cb(1,1) == 2 && cb(2,1) == 2 && model.c == 0)
    error("Condiciones de border invalidas. Neumann-Neumann con c=0");
  endif

  % constantes
  k = model.k;
  c = model.c;
  rho = model.rho;
  c_p = model.c_p;
  G = model.G;
  N = length(xnode)-2;    % nodos interiores
  h = (xnode(end) - xnode(1))/(N+1);

  % matriz base (todo dirichlet, malla uniforme, linda)
  a_i = -(2 + (h^2*c)/k) * ones(N+2,1); % diagonal principal
  b_i = ones(N+1,1);                  % diagonal superior
  c_i = ones(N+1,1);                  % diagonal inferior

  d_i = zeros(N+2,1);
  d_i(2:N+1) = -h^2 * G(xnode(2:N+1)) / k;     % termino independiente (len N)

  M = diag(a_i) + diag(c_i, -1) + diag(b_i, 1);

  % condicion por izquierda
  if (cb(1,1) == 1)
    % dirichlet
    M(1,:) = [1,zeros(1,N+1)];  % phi(1) = q
    d_i(1) = cb(1,2);
  elseif (cb(1,1) == 2)
    % neumann
    % POST DIRICH
  elseif (cb(1,1) == 3)
    % robin
    % POST NEUMANN
  endif

  % condicion por derecha
  if (cb(2,1) == 1)
    % dirichlet
    M(N+2,:) = [zeros(1,N+1),1]; % phi(N+2) = q
    d_i(end) = cb(2,2);
  elseif (cb(2,1) == 2)
    % neumann
    % POST DIRICH
  elseif (cb(1,1) == 3)
    % robin
    % POST NEUMANN
  endif

  % resolver sistema
  T = M \ d_i(:);
endfunction
