o2=[25.48,1.520*(10.)^(-2),-0.7155*(10.)^(-5),1.312*(10.)^(-9),0];
ch4=[19.89,5.024*(10.)^(-2),1.269*(10.)^(-5),-11.01*(10.)^(-9),-72.1*(10.)^3];
co2=[22.26,5.981*(10.)^(-2),-3.501*(10.)^(-5),7.469*(10.)^(-9),-394.0*(10.)^3];
h2o=[32.24,0.1923*(10.)^(-2),1.055*(10.)^(-5),-3.595*(10.)^(-9),-244.5*(10.)^3];
n2=[28.90,-0.1571*(10.)^(-2),0.8081*(10.)^(-5),-2.873*(10.)^(-9),0];

cp(T,x)=T*x[1]+x[2]*(T^2)/2.+x[3]*(T^3)/3.+x[4]*(T^4)/4.;
cp_p(T,x)=x[1]+x[2]*T+x[3]*(T^2)+x[4]*(T^3);

Tref=298;#Temperatura referencia
Ti=Tref;#Temperatura entrada
To=1850;#Tref*np.linspace(2100/298.,2150/298.,101)#Temperatura de salida
#To=linspace(2100/298.,2150/298.,101)

n_ch4_i=1;n_o2_i=2;n_n2_i=7.524;#moles entrada
n_co2_o=1;n_h2o_o=2;n_n2_o=7.524;n_o2_o=1;#moles salida
EA=100;#exceso de aire%100

#1*ch4+2*o2+7.524*n2---->2*h2o+1*co2+7.524*n2-----------
entrada=n_ch4_i*(cp(Ti,ch4)-cp(Tref,ch4)+ch4[5])+(EA/100.)*(n_o2_i*(cp(Ti,o2)-cp(Tref,o2))+n_n2_i*(cp(Ti,n2)-cp(Tref,n2)));#balance entalpias entrada

salida=n_co2_o*(cp(To,co2)-cp(Tref,co2)+co2[5])+n_h2o_o*(cp(To,h2o)-cp(Tref,h2o)+h2o[5])+(EA/100.)*n_n2_o*(cp(To,n2)-cp(Tref,n2))+((EA-100)/100.)*n_o2_o*(cp(To,o2)-cp(Tref,o2));#balance entalpias salida

f(To1)=abs(n_co2_o*(cp(To1,co2)-cp(Tref,co2)+co2[5])+n_h2o_o*(cp(To1,h2o)-cp(Tref,h2o)+h2o[5])+(EA/100.)*n_n2_o*(cp(To1,n2)-cp(Tref,n2))+((EA-100)/100.)*n_o2_o*(cp(To1,o2)-cp(Tref,o2))-entrada);

f_p(To1)=abs(n_co2_o*(cp_p(To1,co2)+co2[5])+n_h2o_o*(cp_p(To1,h2o)+h2o[5])+(EA/100.)*n_n2_o*(cp_p(To1,n2))+((EA-100)/100.)*n_o2_o*(cp_p(To1,o2))-entrada)

arr=map(f,linspace(2000,2100,101));
arr1=linspace(2000,2100,101);

#derivada=(arr[2:end]-arr[1:end-1])./(arr1[2:end]-arr1[1:end-1]))

function derivada(y,x)
arr=linspace(x-1e-2,x+1e-2,3)
return((map(y,arr[2:end])-map(y,arr[1:end-1]) )./(arr[2:end]-arr[1:end-1]))[2]
end

newton(To1)=To1-f(To1)/derivada(f,To1)#newton-raphson

#println(derivada(f,2100.))
#println(minimum(arr))
#println(linspace(2000,2500,10001)[findmin(arr)[2]])

ii=0;i=2000.;
while ii<10
	i=newton(i);
	println(i);#println(ii)
	ii+=1;
end
