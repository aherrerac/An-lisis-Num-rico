%Funcion que muestra las graficas de error vs hi
function tarea2()
  output_precision(16) %Precision double
  [hi,ec,ef,eb]=diff_error(1);
  subplot (1,3,1);
  loglog(hi,ec,hi,ef,hi,eb);
  xlabel ("hi");
  ylabel ("Erel");
  title ("f(x)' / x=1");
  legc = legend ({"Centrado"}, "Adelante","Atras");
  legend (legc, "location", "southwest");
  set (legc, "fontsize", 10);
  [hi,ec,ef,eb]=diff_error(0.1);
  subplot (1,3,2);
  loglog(hi,ec,hi,ef,hi,eb);
  xlabel ("hi");
  ylabel ("Erel");
  title ("f(x)' / x=0.1");
  legf = legend ({"Centrado"}, "Adelante","Atras");
  legend (legf, "location", "southwest");
  set (legf, "fontsize", 10);
  [hi,ec,ef,eb]=diff_error(0.66);
  subplot (1,3,3);
  loglog(hi,ec,hi,ef,hi,eb);
  xlabel ("hi");
  ylabel ("Erel");
  title ("f(x)' / x=0.66");
  legb = legend ({"Centrado"}, "Adelante","Atras");
  legend (legb, "location", "southwest");
  set (legb, "fontsize", 10);
endfunction
%Funcion evaluada
function retval = f(x)
  retval = sin(e^(x^2));
endfunction
%Derivada analitica de la funcion evaluada
function retval = diff_analitic(x)
  retval = cos(e^(x^2))*(e^(x^2))*(2*x);
endfunction
%Diferencia centrada
function retval = diff_centered(x,h)
  retval = (f(x+h)-f(x-h))/(2*h);
endfunction
%Diferencia hacia adelante
function retval = diff_forward(x,h)
  retval = (f(x+h)-f(x))/h;
endfunction
%Diferencia hacia atras
function retval = diff_backward(x,h)
  retval = (f(x)-f(x-h))/h;
endfunction
%Funcion que realiza las aproximaciones y devuelve el error 
function[hi,ec,ef,eb]=diff_error(x)
h=1;           %h inicial
for k = 1:3400 %Muestra de 15 decadas
  a(k,1)=h;
  a(k,2)=diff_centered(x,h);
  a(k,3)=diff_forward(x,h);
  a(k,4)=diff_backward(x,h);
  h*=0.99;     %factor lambda
endfor
hi= a(:,1);
tv = diff_analitic(x); 
ec = abs((tv-a(:,2))/tv);  %Error centrado
ef = abs((tv-a(:,3))/tv);  %Error adelante
eb = abs((tv-a(:,4))/tv);  %Error atras
endfunction



