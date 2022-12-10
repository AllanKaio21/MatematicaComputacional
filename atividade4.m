function atividade4()
  clc;
  x0 = [-0.2; -0.2];
  iter = 100;
  tol = 10e-5;
  metodo = 2;
  if metodo == 1
    nome = 'Metodo do Gradiente';
    a = 0.002;
    [x0, i, todosX, todosY] = metodoGradiente(x0, iter, tol, a);
  elseif metodo == 2
    nome = 'Metodo de Newton';
    a = 0.81;
    [x0, i, todosX, todosY] = metodoNewton(x0, iter, tol, a);
  elseif metodo == 3
    nome = 'Metodo BFGS';
    a = 0.126001;
    [x0, i, todosX, todosY] = metodoBFGS(x0, iter, tol, a);
  endif
  fprintf('\n%s:\nAlfa: %1.5f\nIterações: %d\nResultado: x1 = %1.6f, x2 = %1.6f\nfunção(x1, x2) = %1.6f\n', nome, a, i-1, x0, sistema(x0));
  plotaGraficosConvergencia(i, todosX, todosY);
  plotaGrafico3D(todosX, todosY, i);
endfunction

function [x0, i, todosX, todosY] = metodoGradiente(x0, iter, tol, a)
  todosX(1, 1:2) = x0;
  todosY(1, 1:2) = sistema(x0);
  for i = 2:iter+1
    g = gradiente(x0);
    x1 = x0 - a*g;
    if max(abs(x1 - x0)) <= tol % calcula a tolerancia
      x0 = x1;
      todosX(i + 1, 1:2) = x0;
      todosY(i + 1, 1:2) = sistema(x0);
      break;
    endif
    x0 = x1;
    todosX(i + 1, 1:2) = x0;
    todosY(i + 1, 1:2) = sistema(x0);
  endfor
endfunction

function [x0, i, todosX, todosY] = metodoNewton(x0, iter, tol, a)
  todosX(1, 1:2) = x0;
  todosY(1, 1:2) = sistema(x0);
  for i = 2:iter+1
    h = hessiana(x0);
    g = gradiente(x0);
    x1 = x0 - a*inv(h)*g;
    todosX(i + 1, 1:2) = x1;
    todosY(i + 1, 1:2) = sistema(x1);
    if max(abs(x1 - x0)) <= tol % calcula a tolerancia
      x0 = x1;
      todosX(i + 1, 1:2) = x0;
      todosY(i + 1, 1:2) = sistema(x0);
      break;
    endif
    x0 = x1;
    todosX(i + 1, 1:2) = x0;
    todosY(i + 1, 1:2) = sistema(x0);
  endfor
endfunction

function [x0, i, todosX, todosY] = metodoBFGS(x0, iter, tol, a)
  D = eye(2);
  todosX(1, 1:2) = x0;
  todosY(1, 1:2) = sistema(x0);
  for i = 2:iter+1
    g = gradiente(x0);
    s = -a*D*g;
    y = gradiente(x0+s) - g;
    D =  D + ((s'*y+y'*D*y)*(s*s'))/((s'*y)^2) - (D*y*s' + s*y'*D)/(s'*y);
    d = a*D*g;
    x1 = x0 - d;
    if max(abs(x1 - x0)) <= tol % calcula a tolerancia
      x0 = x1;
      todosX(i + 1, 1:2) = x0;
      todosY(i + 1, 1:2) = sistema(x0);
      break;
    endif
    x0 = x1;
    todosX(i + 1, 1:2) = x0;
    todosY(i + 1, 1:2) = sistema(x0);
  endfor
endfunction

function H = hessiana(x)
  H(1,1) = 2 + 40*pi^2*cos(2*pi*x(1));     % derivada de u em relacao a x1x1
  H(1,2) = 0;                              % derivada de u em relacao a x1x2
  H(2,1) = 0;                              % derivada de v em relacao a x2x1
  H(2,2) = 2 + 40*pi^2*cos(2*pi*x(2));     % derivada de v em relacao a x2x2
end

function G = gradiente(x)
  G(1,1) = 2*x(1) + 20*pi*sin(2*pi*x(1));  % derivada de u em relacao a x1
  G(2,1) = 2*x(2) + 20*pi*sin(2*pi*x(2));  % derivada de u em relacao a x2
endfunction

function y = sistema(x)
    y = 10*(2 - cos(2*pi*x(1)) - cos(2*pi*x(2))) + x(1)^2 + x(2)^2;
endfunction

function plotaGraficosConvergencia(iter, todosX, todosY)
  clf;
  %variaveis
  figure(1);
  subplot(211)
  plot(0:iter, todosX, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado das raizes', 'fontsize', 14);
  title('Grafico de convergencia da variaveis', 'fontsize', 15);
  legend('x1','x2');
  %sistema
  subplot(212)
  plot(0:iter, todosY, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado do sistema', 'fontsize', 14);
  title('Grafico de convergencia do sistema', 'fontsize', 15);
  legend('u(x1,x2)','v(x1,x2)');
end

function plotaGrafico3D(todosX, todosY, iter)
  maiorValor = max(max(todosX));
  menorValor = min(min(todosX));
  intervalo = menorValor:0.01:maiorValor;
  x = intervalo;
  y = x;
  F = zeros(length(x), length(y), 2);
  for xx = 1:numel(x)
    for yy = 1:numel(y)
      F(xx, yy, 1:2) = sistema([x(xx); y(yy)]);
    end
  end
  for cont = 1:iter
    figure(2);
    surf(x, y, F(:, :, 1));
    hold on
    plot3(todosX(1:cont, 2), todosX(1:cont, 1), todosY(1:cont, 1), 'color', 'r', 'linewidth', 2)
    plot3(todosX(cont, 2), todosX(cont, 1), todosY(cont, 1), 'o', 'markersize', 15, 'color', 'r', 'markerfacecolor', [1 1 1], 'linewidth', 2)
    hold off
    title(sprintf('Grafico Animado - iteração %i', cont-1));
    set(gca, 'fontsize', 15);
    xlabel('x1');
    ylabel('x2');
    zlabel('f(x1, x2)');
    grid on;
    view([45 60]);
    pause(0.01);
  endfor
endfunction

