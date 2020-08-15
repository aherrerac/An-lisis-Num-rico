#Intituto Tecnologico de Costa Rica
#Area Academica de Ingenieria en Computadores
#Metodo RungeKutta cuarto orden
#author: Andre Herrera Chacon
#date  : 15/10/2018
function rk4graphs
  EDO = @(x,y) (x*y*y); 
  h=[1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024];   # step values
  # Get x and y values for different h
  [x1,y1] = rungekutta4(EDO,0,1,1,h(1));
  [x2,y2] = rungekutta4(EDO,0,1,1,h(2));
  [x3,y3] = rungekutta4(EDO,0,1,1,h(3));
  [x4,y4] = rungekutta4(EDO,0,1,1,h(4)); 
  [x5,y5] = rungekutta4(EDO,0,1,1,h(5));
  [x6,y6] = rungekutta4(EDO,0,1,1,h(6));
  [x7,y7] = rungekutta4(EDO,0,1,1,h(7));
  [x8,y8] = rungekutta4(EDO,0,1,1,h(8));
  figure(1)
  #Plot h values vs y
  plot(x1,y1,";h=1/8;",x2,y2,";h=1/16;",x3,y3,";h=1/32;",x4,y4,";h=1/64;",x5,y5,";h=1/128;",
  x6,y6,";h=1/256;",x7,y7,";h=1/512;",x8,y8,";h=1/1024;");
  title("y(x) from y'=x*y^2");
  xlabel("x");
  ylabel("y");
  figure(2)
  # Calculate error :  TrueValue - ActualValue /  True Value
  e = [(abs(2-y1(end))/2),(abs(2-y2(end))/2),(abs(2-y3(end))/2),(abs(2-y4(end))/2),(abs(2-y5(end))/2),(abs(2-y6(end))/2),(abs(2-y7(end))/2),(abs(2-y8(end))/2)];
  #Plot h values vs semilog error 
  semilogy(h,e);
  title("Error rungekutta4");
  xlabel("x");
  ylabel("y"); 
endfunction