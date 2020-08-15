#Intituto Tecnologico de Costa Rica
#Area Academica de Ingenieria en Computadores
#Metodo RungeKutta cuarto orden
#author: Andre Herrera Chacon
#date  : 14/10/2018
function[x,y]=rungekutta4(f,xi,xf,y0,h)
  n=(xf-xi)/h;  #Calculate number of steps
  n2 = h/2;     #Save half of h
  y(1) = y0;    #First y equal to y0
  x=xi:h:xf;    #Create x vector with given parameters
  n6 = h/6;     #Save a sixth of h
  for i = 1 : n
    # K: increment based on the slope at the beginning of the interval
    k1= f(x(i),y(i));
    k2= f(x(i) + n2,y(i) + n2 * k1);
    k3= f(x(i) + n2,y(i) + n2 * k2);
    k4= f(x(i)+h,y(i) + h * k3);
    y(i+1) = y(i) + (k1+2*k2+2*k3+k4) * n6; #Get next y value 
  endfor
endfunction



  
    