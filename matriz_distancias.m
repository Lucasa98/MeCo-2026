% matriz de distancias sin bucles
% Dado P (N x 2), calcular D(i,j) = ||p_i - p_j|| con broadcasting (reshape + permute, sumando en la dimension 3)
% calcular con doble for y comparar con tic/toc para N = 2000

N = 2000;

P = 10*rand(N, 2);

% broadcasting
tic();
A = reshape(P, N, 1, 2);
B = permute(A, [2, 1, 3]);
Dif = A - B;
D = sqrt(sum(Dif.^2, 3));
toc()

% doble for

tic();
D = zeros(N, N);
for i = 1 : N
  for j = 1 : N
    D(i,j) = norm(P(i,:) - P(j,:));
  endfor
endfor
toc()

