% particion por condiciones de contorno

N = 5;

K = 10*rand(5,5);

fijos = [2, 4];

libres = setdiff([1,2,3,4,5,1,2], fijos)

K(libres, libres) \ (f(libres) - K(libres, fijos) * u_fijos)

