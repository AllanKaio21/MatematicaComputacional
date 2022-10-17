function atividade2()
  clc;
  x0 = [1; 1];   %x1 e x2 iniciais
  tolerancia = 1e-5;
  [x0, iter, todosX, todosY, f] = NewtonRaphson(x0, tolerancia);
  fprintf('A raiz do sistema f(x1,x2) = {2x1^2 - 4x1x2 - x2^2 = 0\n');
  fprintf('                            {2x2^2 + 10x1 - x1^2 - 4x1x2 - 5 = 0\n');
  fprintf('f(%1.5f, %1.5f) = %1.6f\n',x0(1), x0(2), f(1));
  fprintf('f(%1.5f, %1.5f) = %1.6f\n',x0(1), x0(2), f(2));
  fprintf('Número de iteracoes %d\n', iter);
  todosX = todosX(1:(iter + 1), :);
  todosY = todosY(1:(iter + 1), :);
  plotaGraficosConvergencia(iter, todosX, todosY);
  maiorValor = max(max(todosX));
  menorValor = min(min(todosX));
  intervalo = menorValor:0.1:maiorValor;
  plotaGrafico3D(intervalo, todosX, todosY, iter);
endfunction

function [x0, i, todosX, todosY, y] = NewtonRaphson(x0, tolerancia)
  todosX = zeros(1001, 2);
  todosY = zeros(1001, 2);
  todosX(1, 1:2) = x0;             %valores de x1 e x2
  todosY(1, 1:2) = sistema(x0);  %resultados do sistema
  for i = 1:1000
    y = sistema(x0);
    x1 = x0 - inv(jacobiana(x0))*y;
    todosX(i + 1, 1:2) = x1;
    todosY(i + 1, 1:2) = sistema(x1);
    if max(abs(x1 - x0)) <= tolerancia  %calcula a tolerancia
      x0 = x1;
      break;
    endif
    x0 = x1;
  endfor
endfunction

function todosY = sistema(todosX)
  x1 = todosX(1);
  x2 = todosX(2);

  todosY = zeros(2,1);
  todosY(1) = 2*x1.^2 - 4*x1*x2 - x2.^2;               %u(x1,x2)
  todosY(2) = 2*x2.^2 + 10*x1 - x1.^2 - 4*x1*x2 - 5;   %v(x1,x2)
endfunction

function J = jacobiana(todosX)
  x = todosX(1);   %x1
  y = todosX(2);   %x2
  J = zeros(2,2);

  J(1,1) = 4*x - 4*y;               %derivada de u em relacao a x1
  J(1,2) = - 4*x - 2*y;   %derivada de u em relacao a x2
  J(2,1) = 10 - 4*y - 2*x;          %derivada de v em relacao a x1
  J(2,2) = 4*y - 4*x;               %derivada de v em relacao a x2
end

function plotaGraficosConvergencia(iter, todosX, todosY)
  %variaveis
  figure(1);
  clf;
  plot(0:iter, todosX, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado das raizes', 'fontsize', 14);
  title('Grafico de convergencia da variaveis', 'fontsize', 15);
  legend('x1','x2');

  %sistema
  figure(2);
  clf;
  plot(0:iter, todosY, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado do sistema', 'fontsize', 14);
  title('Grafico de convergencia do sistema', 'fontsize', 15);
  legend('u(x1,x2)','v(x1,x2)');
end

function plotaGrafico3D(x, todosX, todosY, iter)
  y = x;
  F = zeros(length(x), length(y), 2);
  for xx = 1:numel(x)
    for yy = 1:numel(y)
      F(xx, yy, 1:2) = sistema([x(xx); y(yy)]);
    end
  end
  for cont = 1:iter
    figure(3);

    % plot u(x1,x2)
    subplot(211)
    surf(x, y, F(:, :, 1));
    hold on
    plot3(todosX(1:cont, 2), todosX(1:cont, 1), todosY(1:cont, 1), 'color', 'r', 'linewidth', 2)
    plot3(todosX(cont, 2), todosX(cont, 1), todosY(cont, 1), 'o', 'markersize', 15, 'color', 'r', 'markerfacecolor', [1 1 1], 'linewidth', 2)
    hold off;

    set(gca, 'fontsize', 10);
    xlabel('y','fontsize', 16);
    ylabel('x','fontsize', 16);
    zlabel('u(x1,x2)','fontsize', 14);

    % plot u(x1,x2)
    subplot(212)
    surf(x, y, F(:, :, 2));
    hold on;
    plot3(todosX(1:cont, 2), todosX(1:cont, 1), todosY(1:cont, 2), 'color', 'r', 'linewidth', 2)
    plot3(todosX(cont, 2), todosX(cont, 1), todosY(cont, 2), 'o', 'markersize', 15, 'color', 'r', 'markerfacecolor', [1 1 1], 'linewidth', 2)
    hold off;

    set(gca, 'fontsize', 10);
    xlabel('y','fontsize', 16);
    ylabel('x','fontsize', 16);
    zlabel('v(x1,x2)','fontsize', 14);

    grid on;

    pause(0.5);
  endfor
endfunction
