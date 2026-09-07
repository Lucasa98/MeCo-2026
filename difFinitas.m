function [T] = difFinitas(xnode, model, cb, et, T0)
  % DIFERENCIAS FINITAS 1D MALLA UNIFORME

  % xnode: vector de coordenadas nodales
  % model: struct con las constantes del modelo (k, c, rho, c_p, G(x))
  % cb: matriz de 2X3
  %   fila 1: condicion de borde del lado izquierdo
  %   fila 2: condicion de borde del lado derecho
  %   columna 1: tipo de condicion (1=Dirichlet, 2=Neumann, 3=Robin)
  %   columna 2: temperatura, flujo o coef de conveccion h (dependiendo de la condicion de la columan 1)
  %   columna 3: -1 para Dirichlet y Neumann, valor de temperatura externa para condicion Robin.
  % et: indica el esquema temporal: [tipo val1 val2 val3]
  %   - tipo: 0=estacionario, 1=Forward Euler, 2=Backward Euler, 3=Crank-Nicholson
  %   - val1=maxIt, val2=tolError, val3=dt
  %     (si tipo=1 y no se pasa val3, se usa la mitad del dt critico de estabilidad)
  % T0: (opcional) condicion inicial para los esquemas transitorios. Por defecto, ceros.

  % Casos sin solucion (sin Dirichlet en ningun lado)
  % Solo es un problema en el caso estacionario: en el transitorio la condicion
  % inicial fija el nivel de temperatura y el sistema queda bien planteado.
  if (et(1) == 0 && cb(1,1) == 2 && cb(2,1) == 2 && model.c == 0)
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
    M(1,:) = [-(2+h^2*c/k),2, zeros(1,N)];
    d_i(1) = - h^2*G(xnode(1))/k - 2*h*cb(1,2)/k;     % normal exterior. Para normal interior: + 2*h*cb(1,2)/k
  elseif (cb(1,1) == 3)
    % robin
    M(1,:) = [-(2 + (h^2*c/k) + (2*h*cb(1,2)/k)), 2, zeros(1,N)];
    d_i(1) = - h^2*G(xnode(1))/k - 2*h*cb(1,2)*cb(1,3)/k;
  endif

  % condicion por derecha
  if (cb(2,1) == 1)
    % dirichlet
    M(N+2,:) = [zeros(1,N+1),1]; % phi(N+2) = q
    d_i(end) = cb(2,2);
  elseif (cb(2,1) == 2)
    % neumann
    M(end,:) = [zeros(1,N), 2,-(2+h^2*c/k)];
    d_i(end) = - h^2*G(xnode(end))/k - 2*h*cb(2,2)/k; % normal exterior. Para normal interior: + 2*h*cb(2,2)/k
  elseif (cb(2,1) == 3)
    % robin
    M(end,:) = [zeros(1,N), 2, -(2 + (h^2*c/k) + (2*h*cb(2,2)/k))];
    d_i(end) = - h^2*G(xnode(end))/k - 2*h*cb(2,2)*cb(2,3)/k;
  endif

  % CASO ESTACIONARIO: resolver y a pelarse
  if (et(1) == 0)
    T = M \ d_i(:);
    return;
  endif

  % CASOS TRANSITORIOS
  % rho*c_p * dphi/dt = K*phi + Gv       (y en estacionario K*phi + Gv = 0)
  K  =  (k/h^2) * M;
  Gv = -(k/h^2) * d_i(:);

  esDir  = false(N+2,1);
  phiDir = zeros(N+2,1);
  if (cb(1,1) == 1)
    esDir(1)  = true;   phiDir(1)  = cb(1,2);
  endif
  if (cb(2,1) == 1)
    esDir(end) = true;  phiDir(end) = cb(2,2);
  endif
  K(esDir,:) = 0;
  Gv(esDir)  = 0;

  maxIt = et(2);
  tol   = et(3);

  dtcr = 2*rho*c_p / (4*k/h^2 + c);
  if (length(et) >= 4 && et(4) > 0)
    dt = et(4);
  else
    dt = 0.5*dtcr;
  endif

  % Condicion para Forward Euler/explicito
  if (et(1) == 1 && dt > dtcr)
    warning("difFinitas: dt = %g supera el dt critico = %g, Forward Euler va a divergir", dt, dtcr);
  endif

  % condicion inicial
  if (nargin < 5 || isempty(T0))
    T0 = zeros(N+2,1);
  endif
  phi = T0(:);
  phi(esDir) = phiDir(esDir);   % que arranque cumpliendo Dirichlet

  Id = eye(N+2);
  a  = rho*c_p/dt;

  % Armamos la matriz
  if (et(1) == 2)
    % Backward Euler:   (a*I - K) phi^{n+1} = Gv + a*phi^n
    A = a*Id - K;
  elseif (et(1) == 3)
    % Crank-Nicholson:  (a*I - K/2) phi^{n+1} = (a*I + K/2) phi^n + Gv
    A = a*Id - 0.5*K;
    B = a*Id + 0.5*K;
  endif

  % Avanzar
  for n = 1:maxIt
    phiAnt = phi;

    if (et(1) == 1)
      % Forward Euler:  phi^{n+1} = phi^n + (dt/(rho*c_p)) * (K*phi^n + Gv)
      phi = phiAnt + (dt/(rho*c_p)) * (K*phiAnt + Gv);
    elseif (et(1) == 2)
      phi = A \ (Gv + a*phiAnt);
    elseif (et(1) == 3)
      phi = A \ (B*phiAnt + Gv);
    endif

    phi(esDir) = phiDir(esDir);   % se impone dirichlet

    % corte cuando llega a estado estacionario
    if (norm(phi - phiAnt, Inf) < tol)
      break;
    endif
  endfor

  T = phi;
endfunction
