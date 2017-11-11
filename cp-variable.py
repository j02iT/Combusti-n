from scipy import integrate
from scipy import optimize
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


#coeficientes funcion polinomica
#[a,b,c,d,AH]
o2=[25.48,1.520*10**(-2),-0.7155*10**(-5),1.312*10**(-9),0]
ch4=[19.89,5.024*10**(-2),1.269*10**(-5),-11.01*10**(-9),-72.1*10**3]#-72.1*10**3,-74.89*10**3
co2=[22.26,5.981*10**(-2),-3.501*10**(-5),7.469*10**(-9),-394.0*10**3]#-394.0*10**3,-393.6*10**3
h2o=[32.24,0.1923*10**(-2),1.055*10**(-5),-3.595*10**(-9),-244.5*10**3]#-244.5*10**3,-241.88*10**3
n2=[28.90,-0.1571*10**(-2),0.8081*10**(-5),-2.873*10**(-9),0]
#_________________

class cp:
    def __call__(self,T,x):
        cpdt=T*x[0]+x[1]*(T**2)/2.+x[2]*(T**3)/3.
	#cpdt=T*x[0]+x[1]*(T**2)/2.+x[2]*(T**3)/3.+x[3]*(T**4)/4.
        return cpdt

Tref=298.#Temperatura referencia
Ti=Tref#Temperatura entrada
To=1000.#Tref*np.linspace(350/Tref,2550/Tref,1001)#Temperatura de salida

n_ch4_i=1;n_o2_i=2;n_n2_i=7.524;#moles entrada
n_co2_o=1;n_h2o_o=2;n_n2_o=7.524;n_o2_o=0;#moles salida
EA=100;#exceso de aire%100

#Masa molecular
mm_ch4=12.0107+1.0079*4;mm_o2=15.999*2;mm_n2=14.0067*2;
mm_h2o=1.0079*2+15.999;mm_co2=12.0107+15.999*2;
#___________
cp=cp()#llamada a la funcion
#1*ch4+2*o2+7.524*n2---->2*h2o+1*co2+7.524*n2-----------
#Reactivos----
metano=n_ch4_i*(cp(Ti,ch4)-cp(Tref,ch4)+ch4[4]);
oxigeno=(EA/100.)*(n_o2_i*(cp(Ti,o2)-cp(Tref,o2)));
nitrogeno=n_n2_i*(cp(Ti,n2)-cp(Tref,n2));
#Productos----
dioxido_carbono=n_co2_o*(cp(To,co2)-cp(Tref,co2)+co2[4]);
agua=n_h2o_o*(cp(To,h2o)-cp(Tref,h2o)+h2o[4]);
nitrogeno_salida=(EA/100.)*n_n2_o*(cp(To,n2)-cp(Tref,n2));
oxigeno_salida=((EA-100)/100.)*n_o2_o*(cp(To,o2)-cp(Tref,o2));
#_____________________________________________________________


entrada=n_ch4_i*(cp(Ti,ch4)-cp(Tref,ch4)+ch4[4])+(EA/100.)*(n_o2_i*(cp(Ti,o2)-cp(Tref,o2))+n_n2_i*(cp(Ti,n2)-cp(Tref,n2)))#balance entalpias entrada
salida=n_co2_o*(cp(To,co2)-cp(Tref,co2)+co2[4])+n_h2o_o*(cp(To,h2o)-cp(Tref,h2o)+h2o[4])+(EA/100.)*n_n2_o*(cp(To,n2)-cp(Tref,n2))+((EA-100)/100.)*n_o2_o*(cp(To,o2)-cp(Tref,o2))#balance entalpias salida exceso aire

#n_nitro=n_n2_i*2-n_n2_o*2
###
###funcion para encontrar la temperatura de llama,(iterar)
iterar=lambda t_llama:(cp(t_llama,co2)-cp(Tref,co2)+co2[4])+n_h2o_o*(cp(t_llama,h2o)-cp(Tref,h2o)+h2o[4])+n_n2_o*(cp(t_llama,n2)-cp(Tref,n2))-entrada
#optimizacion, con valor inicial, calcular la raiz iterar=0
f=optimize.newton(iterar,2120) 
#temperatura de llama adiabatica, comprobacion
print(f,iterar(f))

"""
t_plot=np.linspace(Ti,f,1001)


fig = plt.figure()
ax = fig.add_subplot(111)

line, = ax.plot(t_plot,abs(iterar(t_plot)), lw=2)

ax.annotate(int(f), xy=(f,0), xytext=(f-100, 4*100000),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
ax.set_xlim(Ti,f+100)
plt.show();
"""

##temp=298.;#kelvin
masa_agua=2000.0;#g
masa_combustible=0.5123#g
contenido_H=3.7#%
Ti=21.669+273.15;Tf=23.947+273.15;incremento=Tf-Ti;
cm_hilo=10-2.6-0.7#10cm
cp_h2o=4.18#kj/kgK
cp_calorimetro=1738.8#J/C
cp_cable=2.3*4.184#J

Q_h2o=incremento*cp_h2o*masa_agua/1000.#kj
Q_bomba=incremento*cp_calorimetro/1000.#kj
Q_cable=cp_cable*cm_hilo/1000.#kj

#print(Q_h2o,Q_bomba,Q_cable)
pcs=(Q_h2o+Q_bomba-Q_cable)*1000/(4.184*masa_combustible)#cal/g

wh=masa_agua*200/18.
pci=pcs/4.18-220*(contenido_H*masa_combustible/100.)

#print(pcs*4.18,pci*4.18)#cal/g

#calculo poder calorifico
#energia_agua=((masa_agua/mm_h2o)*(cp(temp+incremento,h2o)-cp(temp,h2o)))#moles*kj/(kmol)
#energia_cable=0;masa_combustible=1.;
#poder_calorifico=(energia_agua-energia_cable)/(masa_combustible/mm_ch4)#cal/g
#___
#print(energia_agua)
#__________________________

