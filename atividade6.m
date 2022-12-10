function atividade6()
  clc;
  n = 2;
  b = pi/3;
  a = 0;
  tol = 1e-6;
  I = metodoSimpson(a, b, n);
  exata = quad('funcao', a, b);
  fprintf('Para N = %d segmentos a integral tem como resultado %1.6f, sendo o resultado esperado %1.6f.', n, I, exata);
  if abs(exata - I) < tol
    fprintf('\nResultado atingido com %d segmentos.', n);
  else
    fprintf('\nResultado nao atingido com %d segmentos.', n);
  endif
  Nseg = calculaSegmento(a, b, exata, tol);
  fprintf('\nSao necessarios %d segmentos para atingir o resultado.\n', Nseg);
endfunction

function r = metodoSimpson(a, b, n) % 3/8
  s1 = 0;
  s2  = 0;
  s3 = 0;
  r = 0;
  X = linspace(a, b, n+1);
  for i = 2:n
    if mod(i, 3) == 1
      s1 += funcao(X(i));
    elseif mod(i, 3) == 2
      s2 += funcao(X(i));
    else
      s3 += funcao(X(i));
    endif
  endfor
  r = 3*(b - a) * ((funcao(a) + 3*s1 + 3*s2 + 2*s3 + funcao(b))/(8*n));
endfunction

function y = funcao(x)
  y = sin(x);
endfunction

function N = calculaSegmento(a, b, exata, tol)
  r = 0;
  N = 0;
  while abs(exata - r) > tol
    N = N + 1;
    r = metodoSimpson(a, b, N);
  endwhile
endfunction
