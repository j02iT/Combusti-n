#coeficientes funcion polinomica
#[a,b,c,d,AH]
o2=[25.48,1.520*10**(-2),-0.7155*10**(-5),1.312*10**(-9),0]
ch4=[19.89,5.024*10**(-2),1.269*10**(-5),-11.01*10**(-9),-72.1*10**3]#-72.1*10**3,-74.89*10**3
co2=[22.26,5.981*10**(-2),-3.501*10**(-5),7.469*10**(-9),-394.0*10**3]#-394.0*10**3,-393.6*10**3
h2o=[32.24,0.1923*10**(-2),1.055*10**(-5),-3.595*10**(-9),-244.5*10**3]#-244.5*10**3,-241.88*10**3
n2=[28.90,-0.1571*10**(-2),0.8081*10**(-5),-2.873*10**(-9),0]
#________________________________________________________

class cp:
 def __call__(self,T,x):
  #Elegir---------------------------------------
  #cpdt=T*x[0]+x[1]*(T**2)/2.+x[2]*(T**3)/3.#coeficientes              --a,b,c--
  cpdt=T*x[0]+x[1]*(T**2)/2.+x[2]*(T**3)/3.+x[3]*(T**4)/4.#coeficiente --a,b,c,d-- #descomentar
  return cpdt

Tref=298.#Temperatura referencia
Ti=Tref#Temperatura entrada
#
n_ch4_i=1;n_o2_i=2;n_n2_i=7.524;#moles entrada
n_co2_o=1;n_h2o_o=2;n_n2_o=7.524;n_o2_o=0;#moles salida
EA=100;#exceso de aire%100
#___________
cp=cp()#llamada a la funcion
#1*ch4+2*o2+7.524*n2---->2*h2o+1*co2+7.524*n2-----------
entrada=n_ch4_i*(cp(Ti,ch4)-cp(Tref,ch4)+ch4[4])+(EA/100.)*(n_o2_i*(cp(Ti,o2)-cp(Tref,o2))+n_n2_i*(cp(Ti,n2)-cp(Tref,n2)))#balance entalpias entrada

###funcion para encontrar la temperatura de llama,(iterar),salida-entrada
iterar=lambda t_llama:(cp(t_llama,co2)-cp(Tref,co2)+co2[4])+n_h2o_o*(cp(t_llama,h2o)-cp(Tref,h2o)+h2o[4])+n_n2_o*(cp(t_llama,n2)-cp(Tref,n2))-entrada

##derivada numerica
def derivada(y,x):
 arr=[x-1e-2,x-1e-1,x,x+1e-1,x+1e-2];arr2=[y(i) for i in arr];
 dy=[arr2[i]-arr2[i-1] for i in range(1,len(arr))]
 dx=[arr[i]-arr[i-1] for i in range(1,len(arr))]
 y_p=[dy[i]/dx[i] for i in range(0,len(dx))]
 return(sum(y_p)/float(len(y_p)))

##metodo de Newton--------------------- para calcular raices de la ecuacion x(i+1)=x(i)-(f(x)/f'(x))
def newton(x):
 return(x-iterar(x)/derivada(iterar,x))
v_inicial=1500.;tol2=True;#valorinicial

#resolucion--------------------------- iterar 
#y si el resultado es menor que la tolerancia, (solucion) 
while tol2:
 v1=v_inicial;
 v_inicial=newton(v_inicial);
 v2=v_inicial;
 tol=abs(v1-v2)
 if tol > 1e-5:
  tol2=True
 else:
  tol2=False
t_llama=v2

#resultado
print("%.3f Kelvin"%(t_llama))#solucion
