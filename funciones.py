import scipy.integrate as spi
import numpy as np
import pylab as pl
import cmath  # Complex math module
import math
import random
import numpy as np


#####################################################################
#####################################################################
### esclado
#####################################################################
#####################################################################

def escalado_den(miu = 0.0000391, q=1.5,a = 0.162, b = 0.162,dy=0.5,g = 0.1\
              ,r = 1e3,l = 1e4\
             ,m = 0.3,d = 0.5,k = 0.0000002,p = 15800.0,c = 5\
             ,Y0=1e6 ,I0 = 5e4, T0 = 1e5, W0 = 1e4, V0 = 1e5 , t=1):

  miu = miu/t
  g = g/t
  a = (a * Y0) / t

  b = (b * I0) / t
  dy = dy/t
  #N = N/I0
  r = r/(V0 * t)
  l = l/(T0 * t)
  m = m/t
  d = d/t
  k_T = (k * V0)/t
  k_W = (V0 * T0 *k)/(W0 * t)
  k=k_T
  p = (p*W0)/(V0*t)
  c = c/t
  print("t ..",t)


  return {"miu":miu,\
          "q":q,\
          "a":a,\
          "b":b,\
          "dy":dy,\
          "g":g,\
          "r":r,\
          "l":l,\
          "m":m,\
          "d":d,\
          "k_T":k_T,\
          "k_W":k_W,\
          "k":k,\
          "p":p,\
          "c":c}





#####################################################################
#####################################################################
### ecuaciones modelos
#####################################################################
#####################################################################

def Host_0(x,t,a,b,dy,q,g,miu,d,l,k,m,r,c,p):
  dxdt=[a*q*x[1]*(1-x[0])-(g+miu)*x[0]\
        ,b*x[0]*(1-x[1])- dy *x[1]\
        ,l-k*x[4]*x[2]-m*x[2]\
        ,k*x[4]*x[2]-d*x[3]\
        ,p*x[3]+r*x[1]-c*x[4]]
  return dxdt

def Host_0_scale(x,t,a,b,dy,q,g,miu,d,l,k,m,r,c,p , i0, y0):
  dxdt=[a*q*x[1]*((1/i0)-x[0])-(g+miu)*x[0]\
        ,b*x[0]*((1/y0)-x[1])-dy*x[1]\
        ,l-k*x[4]*x[2]-m*x[2]\
        ,k*x[4]*x[2]-d*x[3]\
        ,p*x[3]+r*x[1]-c*x[4]]
  return dxdt

def Host2_0(t,x,a,q,g,miu,b,dy,d,l,k,m,r,c,p,i0, y0):
    I,Y,T,W,V = x
    dxdt=[a*q*Y*(1/i0-I)-(g+miu)*I\
          ,b*I*(1/y0-Y)-dy*Y\
          ,l-k*V*T-m*T\
          ,k*V*T-d*W\
          ,p*W+r*Y-c*V]
    return dxdt
#####################################################################
#####################################################################
### tau leaping 
#####################################################################
#####################################################################


def deng_stock(INPUT ,tau = 0.005,ND = 400,seed_=182, eq_params={},ruleta = False):
    #INPUT = I,Y,T,W,V
    Change=np.zeros((13,5))

    Change[0,:]=([+1 , 0 , 0 , 0 , 0]);
    Change[1,:]=([-1 , 0 , 0 , 0 , 0]);
    Change[2,:]=([-1 , 0 , 0 , 0 , 0]);
    Change[3,:]=([0 , +1 , 0 , 0 , 0]);
    Change[4,:]=([0 , -1 , 0 , 0 , 0]);
    Change[5,:]=([0 , -1 , 0 , 0 , 0]);
    Change[6,:]=([0 , 0 , +1 , 0 , 0]);
    Change[7,:]=([0 , 0 , -1 , +1 , 0]);
    Change[8,:]=([0 , 0 , -1 , 0 , 0]);
    Change[9,:]=([0 , 0 , 0 , -1 , 0]);
    Change[10,:]=([0 , 0 , 0 , 0 , +1]);
    Change[11,:]=([0 , 0 , 0 , 0 , -1]);
    Change[12,:]=([0 , 0 , 0 , 0 , +1]);
    
    def stoc_eqs(INP , Change = Change , eq_params = eq_params):
      V = INP
      Rate=np.zeros((13))
      miu=eq_params["miu"]
      q=eq_params["q"]
      a=eq_params["a"]
      b=eq_params["b"]
      dy=eq_params["dy"]
      g=eq_params["g"]
      r=eq_params["r"]
      l=eq_params["l"]
      m=eq_params["m"]
      d=eq_params["d"]
      k_T=eq_params["k_T"]
      k_W=eq_params["k_W"]
      k=eq_params["k"]
      p=eq_params["p"]
      c=eq_params["c"]  
      i0=eq_params["i0"]
      y0=eq_params["y0"]   
      #reacciones de I
      Propensity=np.array([[a*q*V[1]*(1/i0)
      ,a*q*V[1]*V[0]
      ,(g+miu)*V[0]
      ,a*V[0]*(1/y0)
      ,a*V[0]*V[1]
      ,dy*V[1]
      ,l
      ,k*V[4]*V[2]
      ,m*V[2]
      ,d*V[3]
      ,p*V[3]
      ,c*V[4]
      ,r*V[1]]])

      normal = np.abs(Propensity)*tau

      tau_step = np.random.poisson(normal).transpose()
      #print("tau step",tau_step)
      if ruleta :
        activate_reaction = np.array([np.random.randint(2, size=len(tau_step))]).transpose()
        sum_change = np.sum(tau_step * Change * activate_reaction , axis=0)
      else:
        sum_change = np.sum(tau_step * Change , axis=0)

      #print("cambio",sum_change)
      #print("V",V)
      V=V+sum_change
      V[V< 0 ] = 0

      return V

    def Stoch_Iteration(INPUT):
        lop=0
        T=[0]
        W=[0]
        V=[0]

        I=[0]
        Y=[0]
        print(seed_)
        np.random.seed(seed_)
        random.seed(seed_)
        for lop in Time:
            res = stoc_eqs(INPUT)
            I.append(INPUT[0])
            Y.append(INPUT[1])

            T.append(INPUT[2])
            W.append(INPUT[3])
            V.append(INPUT[4])
            INPUT=res
        return [I,Y,T,W,V]

    Time=np.arange(0.0, ND, tau)
    [IA,YA,TA,WA,VA]=Stoch_Iteration(INPUT)

    tTime=np.array(Time)
    tTA=np.array(TA)[1:,]
    tWA=np.array(WA)[1:,]
    tVA=np.array(VA)[1:,]

    tIA=np.array(IA)[1:,]
    tYA=np.array(YA)[1:,]

    return tTime, tIA , tYA , tTA , tWA , tVA



#####################################################################
#####################################################################
### múltiples simulaciones estocásticas
#####################################################################
#####################################################################

def stoch_simulations(n_sim \
                      , INPUT \
                      , ND = 50 \
                      , tau = 0.005 \
                      , tau_leaping_f = deng_stock \
                      , eq_dif_params  = {} \
                      , scale_params= {"tscale":1,"I0":1,"Y0":1,"T0":1,"W0":1,"V0":1} ):



  tscale,I0,Y0,T0,W0,V0 = scale_params["tscale"]\
                          ,scale_params["I0"]\
                          ,scale_params["Y0"]\
                          ,scale_params["T0"]\
                          ,scale_params["W0"]\
                          ,scale_params["V0"]

  num_iteraciones = n_sim
  
  np.random.seed(182)
  random.seed(182)
  seed_list = [182] + random.sample(range(10000), num_iteraciones -1 )
  for iter_,seed_i in enumerate(seed_list):
    tTime, tIA , tYA ,tTA, tWA, tVA = tau_leaping_f(INPUT = np.array(INPUT) \
                                                ,ND=ND\
                                                ,tau=tau \
                                                ,seed_ = seed_i \
                                                ,eq_params = eq_dif_params)
    if iter_ == 0:
      tIA_fin, tYA_fin , tTA_fin , tWA_fin, tVA_fin = tIA*I0 \
                                                    , tYA*Y0 \
                                                    , tTA*T0 \
                                                    , tWA*W0 \
                                                    , tVA*V0

      sim1 = {"I" : tIA_fin ,"Y" : tYA_fin \
              ,"T" : tTA_fin ,"W" : tWA_fin ,"V" : tVA_fin}

    else:
      tIA_fin = np.vstack([tIA_fin, tIA*I0])
      tYA_fin = np.vstack([tYA_fin , tYA*Y0])
      tTA_fin = np.vstack([tTA_fin , tTA*T0])
      tWA_fin = np.vstack([tWA_fin , tWA*W0])
      tVA_fin = np.vstack([tVA_fin , tVA*W0])

    if iter_ % 10 == 0 :
      print("iteracion ", iter_)

    
  tIA_avg , tYA_avg , tTA_avg , tWA_avg , tVA_avg = tIA_fin.mean(axis=0) \
                                                , tYA_fin.mean(axis=0) \
                                                , tTA_fin.mean(axis=0) \
                                                , tWA_fin.mean(axis=0) \
                                                , tVA_fin.mean(axis=0)

  tIA_std , tYA_std , tTA_std , tWA_std , tVA_std =tIA_fin.std(axis=0) \
                                                , tYA_fin.std(axis=0) \
                                                , tTA_fin.std(axis=0) \
                                                , tWA_fin.std(axis=0) \
                                                , tVA_fin.std(axis=0)

  tIA_std_plus , tYA_std_plus , tTA_std_plus , tWA_std_plus , tVA_std_plus =tIA_avg + tIA_std \
                                                                          , tYA_avg + tYA_std \
                                                                          , tTA_avg + tTA_std \
                                                                          , tWA_avg + tWA_std \
                                                                          , tVA_avg + tVA_std

  tIA_std_minus , tYA_std_minus , tTA_std_minus , tWA_std_minus , tVA_std_minus =tIA_avg - tIA_std \
                                                                              , tYA_avg - tYA_std \
                                                                              , tTA_avg - tTA_std \
                                                                              , tWA_avg - tWA_std \
                                                                              , tVA_avg - tVA_std
  sim_dict = {"time": tTime/tscale \
              , "sim1" : sim1 \
              , "avg" : {"I" : tIA_avg 
                          ,"Y" : tYA_avg \
                          ,"T" : tTA_avg \
                          ,"W" : tWA_avg \
                          ,"V" : tVA_avg} \
              , "std_plus" : {"I" : tIA_std_plus 
                              ,"Y" : tYA_std_plus \
                              ,"T" : tTA_std_plus \
                              ,"W" : tWA_std_plus \
                              ,"V" : tVA_std_plus} \
              , "std_minus" : {"I" : tIA_std_minus 
                              ,"Y" : tYA_std_minus \
                              ,"T" : tTA_std_minus \
                              ,"W" : tWA_std_minus \
                              ,"V" : tVA_std_minus} \
                }

  return sim_dict

#####################################################################
#####################################################################
### Funciones para el cálculo de puntos de equilibrio 
### y números básicos de reproducción
#####################################################################
#####################################################################

def T_tilde( a1, a2):
  # Calculate the discriminant

  discriminant = np.sqrt((a1**2) - (4 * a2))
  T_plus = (a1 + discriminant) / 2
  T_minus = (a1 - discriminant) / 2

  return T_plus, T_minus

def T_tilde_minus( a1, a2):
  discriminant = np.sqrt((a1**2) - (4 * a2))
  T_minus = (a1 - discriminant) * 0.5

  return T_minus


def calculate_a1_a2(g, m, d, p, T0, Rv0):
  # Calculate a1
  a1 = ((g * d) / (p * m)) + T0 * (1 + 1 / Rv0)

  # Calculate a2
  a2 = (T0**2) / Rv0

  return a1,a2

def calculate_a1(g, m, d, p, T0, Rv0):
  # Calculate a1
  a1 = ((g * d)/ (p * m)) + (T0 * (1 + (1 / Rv0)))
  return a1

def calculate_a2(T0, Rv0):

  # Calculate a2
  a2 = (T0**2) / Rv0

  return  a2

def gE(a,E):
  return a*E

def T0_(lambda_may,m):
  return lambda_may/m

def Rv0_fun(T0,k,p,c,m,d):
  return (T0*k*p)/(c*d)

def RvE(T0,T_tilde_minus):
  return T0/T_tilde_minus

def RvE_2(T0, alpha1, alpha2):
  numerator = 2 * T0
  denominator = alpha1 - math.sqrt(alpha1**2 - 4 * alpha2)
  result = numerator / denominator
  return result


def calculate_Tte(T0, Rv):
  return T0 / Rv

def calculate_Vte(c, g, p, Lambd, m, d, Rv):
  return 1/c * (g + ((p * Lambd) / d) * (1 - 1/Rv))

def calculate_Wte(Lambd, m, d, Rv):
  return (Lambd/d) * (1 - 1/Rv)



def Rdo_calc(alpha, q, beta, delta, mu, gamma):
  return (alpha * q * beta) / (delta * (mu + gamma))

def Iti(Rdo, alpha, q, mu, gamma):
  return (Rdo - 1) / (Rdo * (1 + (mu + gamma) / (alpha * q)))


def Yti(Rdo, delta, beta):
  return (Rdo - 1) / (Rdo * (1 + delta / beta))

def Rdo_calc_esc(alpha, q, beta, delta, mu, gamma,i0,y0):
  return (alpha * q * beta ) / (i0*y0 * delta * (mu + gamma))

###
#####################################################################

def T_til_r0_E(Rv0, E , a  = 2e7, m = 0.31, d = 0.1, p = 1e4, T0 = 5000/0.31):
  '''
  (a , m, d, p, T0, Rv0)
  T0, alpha1, alpha2

  '''

  alpha1 = (a * E * (m + d) / (p * m)) + (T0 * (1 + (1 / Rv0)))
  alpha2 = (T0**2) / Rv0
  result = 0.5 * (alpha1 - np.sqrt(alpha1**2 - 4 * alpha2))


  return result

def W_r0_E(Rv0, E , Lambd = 5000 , a  = 2e7, m = 0.31, d = 0.1, p = 1e4, T0 = 5000/0.31):
    alpha1 = (a * E * (m + d) / (p * m)) + (T0 * (1 + 1 / Rv0))
    alpha2 = (T0**2) / Rv0
    T_til = 0.5 * (alpha1 - np.sqrt(alpha1**2 - 4 * alpha2))
    RVE = T0/T_til
    return (Lambd / (m + d)) * (1 - (1/RVE))



