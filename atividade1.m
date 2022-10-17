function atividade1()
  clc;
  iter = 1000;
  tol = 1e-5;
  xl = 0;
  xu = 5;
  intervaloX(:, 1) = xl:0.001:xu;
  xr = inf;
  [xr, iter, XR, FXR] = falsaPos(xl, xu, iter, xr, tol);
  fprintf('\n\nMétodo da Falsa Posição\n\nA raiz da função 16x*sin(x/10)-(37/2) é f(%1.4f) = %1.6f em %d iterações.\n', xr, funcao(xr), iter);
  plotagemGraficos(iter, intervaloX, XR);
  plotaGraficosConvergencia(iter, FXR, XR);
end

function f = funcao(x)
  f = (16*x*sin(x/10)) - (37/2);
end

function [xr, i, XR, FXR] = falsaPos(xl, xu, imax, xr, tol)
  fl = funcao(xl);
  fu = funcao(xu);
  for i = 1:imax
    xr0 = xu - (fu * (xl - xu) / (fl - fu));
    fr = funcao(xr0);
    if fl * fr < 0
      xu = xr0;
      fu = funcao(xu);
    else
      xl = xr0;
      fl = funcao(xl);
    endif
    if abs(xr0-xr) <= tol
      xr = xr0;
      XR(i) = xr;
      FXR(i) = funcao(xr);
      break;
    endif
    xr = xr0;
    FXR(i) = funcao(xr);
    XR(i) = xr;
  endfor
end

function plotagemGraficos(i, intervaloX, XR)
  Y = calculaFuncaoIntervalo(intervaloX);
  for c = 1:i
    figure(1);
    clf;
    plot(intervaloX, Y, 'linewidth', 2, 'color', [0 1 0]);
    hold on;
    plot(XR(c), funcao(XR(c)), 'o', 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    hold off;

    set(gca, 'fontsize', 16);
    grid on;
    xlim([intervaloX(1) intervaloX(end)]);
    ylabel('f(x)');
    xlabel('x');
    title(sprintf('Funcao: f(%1.4f) = %1.5f. Interacao: %i',XR(c), funcao(XR(c)), i));
    legend('f(x) = x^3 + 2*x^2 - 2', 'Posicao estimada da raiz');
    pause(0.5);
  endfor
end

function plotaGraficosConvergencia(qtdInteracoes, FXR, XR)
  %variaveis
  figure(2);
  clf;
  plot(1:qtdInteracoes, FXR, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado das raizes', 'fontsize', 14);
  title('Grafico de convergencia da variaveis', 'fontsize', 15);
  legend('xr');
  %sistema
  figure(3);
  clf;
  plot(1:qtdInteracoes, XR, 'linewidth', 2);
  set(gca, 'fontsize', 16);
  grid on;
  xlabel('Nº de interacoes', 'fontsize', 14);
  ylabel('Valor estimado do sistema', 'fontsize', 14);
  title('Grafico de convergencia do sistema', 'fontsize', 15);
  legend('f(xr)');
end

function Y = calculaFuncaoIntervalo(X)
  [lin, col] = size(X);
  Y = zeros(lin,col);
  for cont = 1:length(X)
    Y(cont) = funcao(X(cont));
  endfor
end
