#Funcion que calcula el valor de la funcion de error de x 
#con c cantidad de cifras significativas
function[v,ev,ea,n]=anpi_erf(x,c)
  if c < 15
    output_precision(15)
    n=1;
    ea=100;
    es=0.5*10^(2-c);          # Umbral Scarborough
    va=anpi_taylor_erf(x,0);  # Valor aproximado anterior
    while abs(ea)>es
      v=anpi_taylor_erf(x,n); 
      ea=((v-va)/v)*100;
      va=v;
      n++;
    endwhile
    ev=(1-(v/erf(x)));
  else
    display("El numero de cifras significativas debe ser menor que 15");
  endif
endfunction

#Funcion auxiliar que calcula la serie de taylor de la funcion de error
function[s]=anpi_taylor_erf(x,k)
s=0;
  for i = 0 : k
    s+=(((-1)^i)/factorial(i))*((x^(2*i+1))/(2*i+1)); 
  end
s=s*(2/sqrt(pi));
endfunction
#Funcion que devuelve los valores de la cuadratica utilizando la formula 
#alternativa y la tradicional, asi como precision simple y doble
function[x1ts,x2ts,x1as,x2as,x1td,x2td,x1ad,x2ad]=anpi_cuadratica(a,b,c)
    [x1ts,x2ts,x1as,x2as]=anpi_cuadratica_simple(a,b,c);
    display("Precision simple");
    display("Tradicional");
    disp(["x1 =",disp(x1ts)]);
    disp(["x2 =",disp(x2ts)]);
    disp("Alternativa");
    disp(["x1 =",disp(x1as)]);
    disp(["x2 =",disp(x2as)]);
    [x1td,x2td,x1ad,x2ad]=anpi_cuadratica_doble(a,b,c);
    disp("Precision doble");
    display("Tradicional");
    disp(["x1 =",disp(x1td)]);
    disp(["x2 =",disp(x2td)]);
    disp("Alternativa");
    disp(["x1 =",disp(x1ad)]);
    disp(["x2 =",disp(x2ad)]);
endfunction

#Funcion auxiliar que calcula la cuadratica con precision simple
function[x1t,x2t,x1a,x2a]=anpi_cuadratica_simple(a,b,c)
  output_precision(7)
  x1t=single(((-single(b))+sqrt((single(b))^2-4*single(a)*single(c)))/(2*single(a)));
  x2t=single(((-single(b))-sqrt((single(b))^2-4*single(a)*single(c)))/(2*single(a)));
  x1a=(-2*single(c))/(single(b)+sqrt((single(b))^2-4*single(a)*single(c)));
  x2a=single((-2*single(c))/(single(b)-sqrt(single(b)^2-4*single(a)*single(c))));
endfunction
#Funcion auxiliar que calcula la cuadratica con precision doble
function[x1t,x2t,x1a,x2a]=anpi_cuadratica_doble(a,b,c)
  output_precision(15)
  x1t=((-b)+sqrt((b)^2-4*a*c))/(2*a);
  x2t=((-b)-sqrt((b)^2-4*a*c))/(2*a);
  x1a=(-2*c)/(b+sqrt((b)^2-4*a*c));
  x2a=(-2*c)/(b-sqrt((b)^2-4*a*c));
endfunction

