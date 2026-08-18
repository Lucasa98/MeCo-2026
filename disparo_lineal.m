function [x,y] = disparo_lineal(f, inter, yc, L)
  % disparo lineal

  % sistema
  F = @(x,u) [
    u(3);
    u(4);
    f(x)(1)*u(3) + f(x)(2)*u(1) + f(x)(3);
    f(x)(1)*u(4) + f(x)(2)*u(2)
  ];

  % condiciones iniciales y1(a) = alpha, y2(a) = 0, y1'(a) = 0, y2'(a) = 1
  u0 = [
    yc(1);
    0;
    0;
    1
  ];

  puntos = linspace(inter(1), inter(2), L);

  [x, U] = ode45(F, puntos, u0);

  lambda = (yc(2) - U(end, 1))/U(end,2);

  y = U(:,1) + lambda*U(:,2);
endfunction
