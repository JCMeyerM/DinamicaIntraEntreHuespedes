ruta_escritura = "/content/drive/MyDrive/tesis/dengue/simulaciones/"




### parámetros para las simulaciones

#############################################
#############################################
### parámetros simulación determinística
#############################################
#############################################


param_det = {"miu" : 0.0000391
            ,"q":1.5
            ,"a" : 0.162
            ,"b" : 0.162
            ,"dy" : 0.5
            ,"g" : 0.05
            ,"r" : 1e3
            ,"l" : 1e4
            ,"m" : 0.3
            ,"d" : 0.5
            ,"k" : 0.0000002
            ,"p" : 15800.0
            ,"c" : 5}


## escalamiento sol determinista
param_scale_det_0 ={"Y0" : 1 \
                ,"I0" : 1 \
                ,"W0" : 1e4 \
                ,"T0" : 1e4 \
                ,"V0" : 1e5 \
                ,"tscale" : 0.5}


## escalamiento para comprar con la sim estocastica
param_scale_det ={"Y0" : 100 \
                ,"I0" : 100 \
                ,"W0" : 1e4 \
                ,"T0" : 1e4 \
                ,"V0" : 1e5 \
                ,"tscale" : 0.5}


# set the initial conditions
ci_det = {"i_ini" : 0
          ,"y_ini" : 0.0001
          ,"T_ini" : 100000
          ,"W_ini" : 0
          ,"V_ini" : 0.1}

tend_det=400 #limite de tiempo
tstep_det= 0.0005 #paso de tiempo

#x0=[i_ini,y_ini,T_ini,W_ini,V_ini]
x0_det=[ci_det["i_ini"],ci_det["y_ini"],ci_det["T_ini"],ci_det["W_ini"],ci_det["V_ini"]]

# valor de mosquitos infectados Y para cálculos de punto s de equilibrio
Y_graf = 0.01


#############################################
#############################################
### parámetros simulación estocástica
#############################################
#############################################

### estocástica
#scale estocasticp
param_stoch = {"Y0" : 0.01
,"I0" : 0.01
,"W0" : 1e2 
,"T0" : 1e2
,"V0" : 1e2
,"tscale" : 0.5}



ci_stoch = {"i_ini" : 0.1\
          ,"y_ini" : 0.1\
          ,"T_ini" : 1e4/0.3\
          ,"W_ini" : 10\
          ,"V_ini" : 10}




stoch_days = 400
tau_step = 0.0005
stoch_seed= 182


### varias simulaciones
num_simulaciones = 50


#############################################
#############################################
### parámetros gráficas R0
#############################################
#############################################

param_r0_2 = {"lambda_may" : 5000 \
            ,"m" : 0.31 \
            ,"d" : 0.1 \
            ,"p" : 1e4 \
            ,"c" : 20 \
            ,"a" : 2e7 \
            ,"k" : 1e-7 } 

param_r0_09 = {"lambda_may" : 5000 \
            ,"m" : 0.31 \
            ,"d" : 0.1 \
            ,"p" : 1e4 \
            ,"c" : 20 \
            ,"a" : 2e7  } 


#############################################
#############################################
### parámetros gráficas 3D
#############################################
#############################################
param_graf3D={"lambda_may" : 5000\
        ,"m" : 0.31\
        ,"d" : 0.1\
        ,"p" : 1e4\
        ,"c" : 20\
        ,"a" : 2e7}


#############################################
#############################################
### parámetros gráficas R0
#############################################
#############################################



