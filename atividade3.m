function atividade3()
  clc;
  fprintf('\nFatoração LU - Pivoteamento Parcial\n');
  fprintf('\nSistema: f(x, y, z) = {-x + y - 3z = -4\n');
  fprintf('                      {3x - 2y + 8z = 14\n');
  fprintf('                      {2x - 2y + 5z = 7\n\n');
  [A, b] = Msistema();
  [L, U] = fatoracaoLU(A);
  d = subProgressiva(L,b);
  x = subRegressiva(U, d);
  fprintf('\nSolução: x = %.2f; y = %.2f; z = %.2f;\n\n', x(:));
endfunction

function [L, U] = fatoracaoLU(U)
  n = size(U);
  L = eye(n(1),n(1)); % Matriz Identidade
  for j = 1:n(1)-1
      for i = j+1:n(1)
          fprintf("\ni: %d j; %d",i,j);
          fator  = U(i,j) / U(j,j);
          L(i,j) = fator;
          U(i,j:n(1)) = U(i,j:n(1)) - (fator*U(j,j:n(1)));
      endfor
  endfor
  L
  U
endfunction

function A = pivoteamento(A)
  n = size(A);
  B = zeros(n(1), n(2));
  for i = 1:n(1)-1
    [~, pos] = max(abs(A(:, i)));
    if(pos != i)
      B(i, :) = A(i,:);
      A(i,:) = A(pos, :);
      A(pos, :) = B(i, :);
    endif
  endfor
endfunction

function d = subProgressiva(L, b)
  Lb = [L b];
  n = size(Lb);
  d = zeros(n(1), 1);
  for i = 1:n(1)
    d(i) = (Lb(i, n(1)+1) - sum(Lb(i,1:i-1)))/Lb(i,i);
    Lb(:, i) = Lb(:, i) * d(i);
  endfor
endfunction

function x = subRegressiva(U, d)
  n = size(U);
  x = zeros(n(1), 1);
  Ud = [U d];
  j = n(1)+1;
  for i = 1:n(1)
    if(j-i==j-1)
      x(j-i) = Ud(j-i, j) / Ud(j-i, j-1);
      Ud(1:j-i, j-i) = Ud(1:j-i, j-i) * x(j-i);
    else
      x(j-i) = (Ud(j-i, j) - sum(Ud(j-i, j-i+1: j-1))) / Ud(j-i, j-i);
      Ud(1:j-i, j-i) = Ud(1:j-i, j-i) * x(j-i);
    endif
  endfor
endfunction

function [A, b] = Msistema()
  A = [-1, 1, -3, -4;
        3, -2, 8, 14;
        2, -2, 5, 7];

  fprintf('Matriz do sistema: \n\n');
  n = size(A,1);
  for i = 1:n
    fprintf('| ');
    fprintf('%.3f  ', A(i,:));
    fprintf('|');
    fprintf('\n');
  endfor

  A = pivoteamento(A);
  b = A(:,size(A,1)+1);
  A = A(:,1:size(A,1));
  # x=4, z=1, y=3

endfunction
