from parameter_values import *

#alpha = 10
#Kd = 1
#n = 2
#delta = 1


def YES(x, params = (alpha, Kd, n)):    
    alpha, Kd, n = params

    frac = (x/Kd)**n
    dy_dt = alpha * (frac)/(1 + frac)
    
    return dy_dt

def NOT(x, params = (alpha, Kd, n)):    
    alpha, Kd, n = params

    frac = (x/Kd)**n
    dy_dt = alpha * (1)/(1 + frac)
    
    return dy_dt

def AND(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (frac1*frac2)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt


def AND00(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (1)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt


def AND01(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (frac2)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt

def AND10(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (frac1)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt

def AND11(in1, in2, params = (alpha, Kd, n)):        
    
    return AND(in1, in2, params)




def OR(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (frac1 + frac2 + frac1*frac2)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt

def NOR(in1, in2, params = (alpha, Kd, n)):        
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    dy_dt = alpha * (1)/(1 + frac1 + frac2 + frac1*frac2)
    
    return dy_dt


def EQU(in1, in2, params = (alpha, Kd, n)):
    alpha, Kd, n = params

    frac1 = (in1/Kd)**n
    frac2 = (in2/Kd)**n
    
    gate1 = (frac1*frac2)/(1 + frac1 + frac2 + frac1*frac2) # and
    gate2 = (1)/(1 + frac1 + frac2 + frac1*frac2) # nor
    dy_dt = alpha * (gate1 + gate2)

    #x1 = AND(in1, in2, params=params)
    #x2 = NOR(in1, in2, params=params)
    #dy_dt=OR(x1, x2, params=params) # wrong - uses alpha**2 instead of alpha

    return dy_dt


def degrade(x, delta=delta):
    dx_dt = - x * delta
    
    return dx_dt