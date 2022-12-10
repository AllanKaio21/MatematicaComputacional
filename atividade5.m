function atividade5()
  clc;
  ordem = 5; # ordem do polinomio
  x = [0.1, 0.2, 0.4, 0.6, 0.9, 1.3, 1.5, 1.7, 1.8];
  y = [0.75, 1.25, 1.45, 1.25, 0.85, 0.55, 0.35, 0.28, 0.18];
  xx = 1.4;
  a = metodo(x, y, ordem);
  yy = funcao(a, xx);
  [r, r2] =  calculaR(y, x, a);

  fprintf("\nAjuste de curvas e interpolação\n");
  fprintf("\nO polinomio fitado e f(x) = ");
  n = size(a, 1);
  for i = 1:n
    fprintf("%1.6f*x^%d",a(n-i+1), n-i);
    if i != n
      fprintf(" + ");
    endif
  endfor
  fprintf("\nsendo f(%1.6f) = %1.6f tendo r^2 = %1.6f e r = %1.6f\n", xx, funcao(a, xx), r2, r);

  plotaGrafico(x, y, yy, xx, a);
endfunction

function a = metodo(x, y, ordem)
  n = size(x, 2);
  for i = 1:ordem+1
    for j = 1:i
      k = i + j - 2;
      soma = 0;
      for w = 1:n
        soma = soma + x(w)^k;
      endfor
      A(i, j) = soma;
      A(j, i) = soma;
    endfor
    soma = 0;
    for w = 1:n
      soma = soma + y(w)*x(w)^(i-1);
    endfor
    A(i, ordem+2) = soma;
  endfor
  a = A(:,1:ordem+1);
  b = A(:, ordem+2);
  a = a\b;
endfunction

function plotaGrafico(x, y, yy, xx, a)
  intervalo = x(1):0.01:x(end);
  todosY = vetFuncao(intervalo, a);
  figure(1);
  clf;
  plot(intervalo, todosY, 'linewidth', 2);
  hold on
  plot(x, y, 'o', 'markerfacecolor', [0 1 0]);
  plot(xx, yy, 'o', 'markerfacecolor', [1 0 0]);
  hold off
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('x', 'fontsize', 14);
  ylabel('f(x)', 'fontsize', 14);
  title('Grafico de Ajuste de curvas e interpolação', 'fontsize', 15);
  legend('Polinomio','Amostras', 'Valor Estimado');
end

function [r, r2] = calculaR(y, x, a)
  st = sum((y - mean(y)).^2);
  vf = vetFuncao(x, a);
  sr = sum((y - vf').^2);
  r2 = (st-sr)/st;
  r = sqrt(r2);
endfunction

function Y = vetFuncao(x, a)
  n = size(x, 2);
  Y = zeros(n, 1);

  for i = 1:n
    Y(i) = funcao(a, x(i));
  endfor
endfunction

function y = funcao(a, x)
  y = 0;
  n = size(a, 1);
  for i = 0:n-1
      y = y + (a(i+1)*x^(i));
  endfor
endfunction
