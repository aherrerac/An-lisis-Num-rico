#Intituto Tecnologico de Costa Rica
#Area Academica de Ingenieria en Computadores
#Metodo RungeKutta cuarto orden
#author: Kevin Acuna Mena 
#date  : 16/10/2018
function EDOMethodComparison
  EDO = @(x,y) (100 - y);
  xi = 0;
  xf = 200;
  y0 = 5;
  h = xf/1000; %Value for h assuring at least 1000 steps.
  
  %RK4 Time and Step Quantity Disposal.
  disp('Runge Kutta:');
  tic
  [xRK4, yRK4] = rungekutta4(EDO,xi,xf,y0,h);
  toc
  msg = ['Puntos utilizados: ', num2str(length(xRK4))]; 
  disp(msg);
  
  %ode45 Time and Step Quantity Disposal.
  disp('ode45:');
  tic
  [xODE45,yODE45] = ode45(EDO, [xi,xf], y0);
  toc
  msg = ['Puntos utilizados: ', num2str(length(xODE45))]; 
  disp(msg);
  
  %ode23 Time and Step Quantity Disposal.
  disp('ode23:');
  tic
  [xODE23,yODE23] = ode23(EDO, [xi,xf], y0);
  toc
  msg = ['Puntos utilizados: ', num2str(length(xODE23))]; 
  disp(msg);
  
  %RK4 Column Vector Construction for 100<x<200 and 99.8<y<100.2.
  x1 = 0;
  y1 = 0;
  counter = 1;
  for i = 1:length(xRK4)
    if((xRK4(i) >= 100)&&(yRK4(i) >= 99.8))
      x1(counter,1) = xRK4(i);
      y1(counter,1) = yRK4(i);
      ++counter;
    endif    
  endfor
  
  %ode45 Column Vector Construction for 100<x<200 and 99.8<y<100.2.
  x2 = 0;
  y2 = 0;
  counter = 1;
  for i = 1:length(xODE45)
    if((xODE45(i) >= 100)&&(yODE45(i) >= 99.8))
      x2(counter,1) = xODE45(i);
      y2(counter,1) = yODE45(i);
      ++counter;
    endif
  endfor
  
  %ode23 Column Vector Construction for 100<x<200 and 99.8<y<100.2.
  x3 = 0;
  y3 = 0;
  counter = 1;
  for i = 1:length(xODE23)
    if((xODE23(i) >= 100)&&(yODE23(i) >= 99.8))
      x3(counter,1) = xODE23(i);
      y3(counter,1) = yODE23(i);
      ++counter;
    endif
  endfor
  
  %Graph Disposal.
  plot([100;x1;200],[99.8;y1;100.2],";RK4;",[100;x2;200],[99.8;y2;100.2],";ode45;",[100;x3;200],[99.8;y3;100.2],";ode23;");
  title("Method Comparison for y' = 100 - y");
  xlabel("x");
  ylabel("y"); 
endfunction