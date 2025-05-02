# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 18:54:59 2024

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 21:48:06 2024

@author: user
"""

import sympy as sy
import numpy as np
import pylab as plt
import pandas as pd
import math
import matplotlib.animation as animation
import sys

"Nomes de perfis que tem no arquivo"
# AeroJampa 2021
# Selig 1223
# FX 63-0137 13.7%
# AeroJampa 2022
# Selig 1210
# Selig 1223 RTL
# EPPLER 420
# EPPLER 423
# NACA 4412
# NACA 4415
# FX 73-CL3-152
# NACA 0015
# NACA 63015a
# NACA 0009
# NACA 0012
# NACA 64012
# BELL 540


"inputs Físicos"

Temp = 20           #[°C] Temperatura que você quer analisar a decolagem
# g   = 9.7807        #[m/s²] Aceleração gravitacional  em João pessoa
g   = 9.7896        #[m/s²] Aceleração gravitacional  em São José dos campos


"inputs Aerodinamica"

h_asa   = 0.28       #[m] altura da asa em relação ao chão

i_w = 5          #[graus] ângulo de incidência da asa
torção_asa = 0   #[graus] ângulo de torção da asa
diedro_asa = 0   #[graus] ângulo de diedroda asa
alpha_s = 20     #[graus] ângulo de estol da asa


Sep     = 0.2                  #[m²] área lateral de 1 endplate, caso não tenha endplate, coloque esse valor = 0
h_end   = 0.4                  #[m] altura maxima de 1 endplate, caso não tenha endplate, coloque esse valor = 0
Swetend = 0.4                  #[m²] área molhada dos endplates(eixo x)
Ced     = 0.5                  #[m] corda/raio do endplate
t_end   = 0.4                  #Espessura maxima do endplate (geometria)
Cd0ed   = 0.9                  #Cd0 do endplate (geometria)

Cd0_f   = 1.2                  #Cd0 da fuselagem (geometria)

perfil_asa = "Selig 1223 RTL"
perfil_EH  = "NACA 0015"
perfil_EV  = "NACA 0015"


lista_Cr_asa = [0.55,0.55]   #corda raiz da primeira seção, segunda seção, etc...
lista_Ct_asa = [0.55,0.3]    #corda da ponta da primeira seção, segunda seção, etc...
lista_b_asa  = [0.55,1.05]   #posição em b(eixo y) do fim da primeira seção, segunda seção, etc...

#lembrar que os valores de b na lista é divido por 2

            
#Fuselagem

lf=0.7873        #Comprimento da fuselagem eixo x
df=0.133         #largura da fuselagem eixo y 
         

"inputs Controle"

S_EH = 0.223     #[m²] Área do EH
S_EV = 0.043     #[m²] Área do EV
i_EH = -1        #[graus] ângulo de incidência do EH

lista_Cr_EH = [0.275]   #corda raiz da primeira seção, segunda seção, etc...
lista_Ct_EH = [0.275]   #corda da ponta da primeira seção, segunda seção, etc...
lista_b_EH  = [0.405]   #posição em b(eixo y) do fim da primeira seção, segunda seção, etc...
            
#lembrar que os valores de b na lista é divido por 2

lista_Cr_EV = [0.175]   #corda raiz da primeira seção, segunda seção, etc...
lista_Ct_EV = [0.175]   #corda da ponta da primeira seção, segunda seção, etc...
lista_b_EV  = [0.1388]  #posição em b(eixo y) do fim da primeira seção, segunda seção, etc...

#lembrar que os valores de b na lista é divido por 2

"inputs Performace"

"Lembre de mudar a enquação de tração do motor la embaixo"

atrito = 0.0658  #[adimensional] coeficiente de atrito entre as rodas e o chão
pista_max = 55   #[m] Posicao do obstáculo em x
h_obs = 0.7      #[m] altura do obstaculo
x_ativa_profundor = pista_max* 0.5 #[m] afeta somente analise ideal, momento que o profundor é acionado

"inputs do Código(afetam a acuracidade dos resultados e tempo de processamento"

comecar_MTOW     = 5
MTOW_final       = 25
incremento_MTOW  = 0.1
incremento_tempo = 0.01

"Analise ""ideal"" considera que na metade da pista, o profundor é acionado e levanta o nariz do avião até-"
"- o ângulo de estol da asa, alcançando CLmax"

"Analise ""conservativa"" considera o avião sem nenhuma alteração de controle, onde somente quando a asa -"
"gera sustentação suficiente para levantar o avião ele muda seu ângulo, chegando no CLmax"

analise = "ideal" #"ideal" ou "conservativa"

"Aqui começa___________________________________________________________________________________"

Vx = sy.symbols('Vx')
Vy = sy.symbols('Vy')
V_t = sy.symbols("V_t")
m = sy.symbols('m')
alpha = sy.symbols('alpha')



"Equação de tração do motor"

T = -0.05816*(V_t**2) + 0.4061*V_t + 41.47  #equação de tração do motor em função da velocidade

"Arquivos dos dados de perfil"
"O arquivo tem valores como cd0, Clalpha, Cl0, Clmax, a_s etc em função do número de Reynolds dividido por 10000"

arquivos_perfil = pd.read_excel("Perfil_df.xlsx")


rho = -3.52607427278736e-8*Temp**3 + 1.60215687088211e-5*Temp**2 - 0.00471578515074178*Temp + 1.2925545912127
print(f"rho em função da temperatura = {rho}")
rho = 1.1116
print(f"rho usado = {rho}")
vis = 5.63034007726774e-13*Temp**3 + 1.15655321938951e-11*Temp**2 + 9.06376414877547e-8*Temp + 1.32700903275642e-5


if analise != "conservativa" and analise != "ideal":
    print("Esolha entre analise conservativa ou ideal!")
    sys.exit()


def integral_sem_integral(vx,vy):
    
    area = 0

    for i in range(len(vx)-1):
        
        h_tri = abs(vy[i] - vy[i+1])
        
        base  =  vx[i+1] - vx[i]
        
        # h_ret = min([(vy[i]),vy[i+1]])
        h_ret = min([abs(vy[i]),abs(vy[i+1])])

        a_tri =  (base * h_tri)/2
        
        a_ret = base * h_ret
        
        # ou
        
        #a_trap = (min([vy[i],vy[i+1]]) +(min([vy[i],vy[i+1]])+abs(vy[i] - vy[i+1])))*(vx[i+1] - vx[i])/2

        area = area + a_tri + a_ret
    
    return(area)
    
    

def calculo_MAC(lista_Cr,
                lista_Ct,
                lista_b):
    
    'calculo da corda média'
    
    lista_S = []
    lista_C = []
    
    for i in range(len(lista_b)):
        
        TR_i = lista_Ct[i]/lista_Cr[i]
        cmac=(2/3)*lista_Cr[i]*((1+TR_i+(TR_i**2))/(1+TR_i))
        
        if i == 0:
            S_i = lista_b[i]*cmac
        if i > 0:
            S_i = (lista_b[i]-lista_b[i-1])*cmac
            
        
        S_i = 2*S_i
        
        lista_C.append(cmac)
        lista_S.append(S_i)
    
    soma = sum(lista_S)
    soma_cima = 0
    for i in range(len(lista_C)):
        
        soma_cima = soma_cima + (lista_C[i]*lista_S[i])
        
    
    MAC = soma_cima/soma
    
    return(MAC)

def LLT(lista_b:list,
        lista_Cr:list,
        lista_Ct:list,
        Clalpha:float,
        torção:float,
        incidencia:float,
        alpha_0:float,
        diedro:float,
        N:int):
    
    "Essa função transforma os dados de perfil de Cl e transforma em valores para asa 3D CL"
    
    twist = torção
    Cla = Clalpha
    i_w = incidencia
    
    '__________________________________________________________________________'
    
    'Calculo das linhas que representam a asa'
    
    x=sy.symbols('x')
    
    segmentos=[]
    
    for i in range(len(lista_b)):
        if i == 0:
            linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-0))*(x-0))
        if i > 0:     
            linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-lista_b[i-1] + 0.0001))*(x-lista_b[i-1]))
        
        segmentos.append(linha)
    
    listax = np.linspace(0,lista_b[len(lista_b)-1],N) #divido N vezes
    listay = []
    listacont = np.arange(1,len(lista_b),1)
    
    for x in listax:
        
        if x < lista_b[0]:
            y = eval(str(segmentos[0]))
        
        if x > lista_b[0]:
            for i in listacont:
                if x > lista_b[i-1] and x < lista_b[i]:
                    y = eval(str(segmentos[i]))
                    
        listay.append(y)
    
    '___________________________________________________________________________'
    
    'calculo da corda média'
    
    lista_S = []
    lista_C = []
    
    for i in range(len(lista_b)):
        
        TR_i = lista_Ct[i]/lista_Cr[i]
        cmac=(2/3)*lista_Cr[i]*((1+TR_i+(TR_i**2))/(1+TR_i))
        
        if i == 0:
            S_i = lista_b[i]*cmac
        if i > 0:
            S_i = (lista_b[i]-lista_b[i-1])*cmac
            
        
        S_i = 2*S_i
        
        lista_C.append(cmac)
        lista_S.append(S_i)
    
    soma = sum(lista_S)
    soma_cima = 0
    for i in range(len(lista_C)):
        
        soma_cima = soma_cima + (lista_C[i]*lista_S[i])
        
    
    MAC = soma_cima/soma
    
    b = lista_b[len(lista_b)-1]*2
    
    Cla = Cla * 180/np.pi
    
    alpha_twist = twist         
    a_2d = Cla
    
    theta = np.linspace((np.pi / (2 * N)), (np.pi / 2), N, endpoint=True)
    
    alpha = np.linspace(i_w + alpha_twist, i_w, N) 

    listay = listay[::-1]
    
    c = np.array(listay)
    
    mu = c * a_2d / (4 * b)
    
    LHS = mu * (np.array(alpha) - alpha_0) *np.pi/180
    
    RHS = []
    for i in range(1, 2 * N + 1, 2):
        RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))  # .reshape(1,self.N)
        # print(RHS_iter,"RHS_iter shape")
        RHS.append(RHS_iter)
        
    
    test = np.asarray(RHS)
    x = np.transpose(test)
    inv_RHS = np.linalg.inv(x)
        
    A = np.matmul(inv_RHS, LHS) #An
        
    AR = b/MAC         #Aspect ratio
    
    CL_wing = (np.pi * AR * A[0]) 
    
    CL_wing = CL_wing * np.cos(diedro*np.pi/180)
    
    listay = listay[::-1]  
    
    return(round(CL_wing,5),listax,listay)

def L_e_D(arquivos_perfil,
            #asa
            perfil_asa:str,
            lista_b_asa:list,
            lista_Cr_asa:list,
            lista_Ct_asa:list,
            torção_asa:list,
            diedro_asa:list,
            i_w:float,
            h:float,
            alpha:float,
            
            #EH
            perfil_EH:str,
            lista_b_EH:list,
            lista_Cr_EH:list,
            lista_Ct_EH:list,
            
            #EV
            perfil_EV:str,
            lista_b_EV:list,
            lista_Cr_EV:list,
            lista_Ct_EV:list,
            
            #Fuselagem
            lf:float,
            df:float,
            
            #Constantes
            V:float,
            rho:float,
            vis:float,
            ):
    
    if V < 0 :
        return(0,0)
    
    '________________________________________________________________________________________'

    
    MAC_asa = calculo_MAC(lista_Cr_asa, lista_Ct_asa, lista_b_asa)
    MAC_EH =  calculo_MAC(lista_Cr_EH, lista_Ct_EH, lista_b_EH)
    MAC_EV =  calculo_MAC(lista_Cr_EV, lista_Ct_EV, lista_b_EV)
    
    Re   = V*MAC_asa/vis
    Reht = V*MAC_EH/vis
    Revt = V*MAC_EV/vis
    Ref  = V*lf/vis
    Rend = V*Ced/vis

    
    nomes    = list(arquivos_perfil["Nome"])
    clalphas = list(arquivos_perfil["Clα vs Re"])
    cd0s     = list(arquivos_perfil["Cd0 vs Re"])
    alphas_0 = list(arquivos_perfil["α0 vs Re"])
    a_sl     = list(arquivos_perfil["αs vs Re"])
    t_cs     = list(arquivos_perfil["t/c"])
    
    
    x = Re/10000
    
    if x < 10:
        x = 10
    
    index_w   = nomes.index(perfil_asa)
    Clalpha_w = eval(clalphas[index_w])
    Cd0_w     = abs(eval(cd0s[index_w]))
    alpha0_w  = eval(alphas_0[index_w])
    a_s       = eval(a_sl[index_w])
    t_c_w     = t_cs[index_w]
    
    x = Reht/10000
    
    if x < 10:
        x = 10
    
    index_ht   = nomes.index(perfil_EH)
    Clalpha_EH = eval(clalphas[index_ht])
    Cd0_EH     = abs(eval(cd0s[index_ht]))
    alpha0_EH  = eval(alphas_0[index_ht])
    a_s_EH     = eval(a_sl[index_ht])
    t_c_EH     = t_cs[index_ht]
    
    x = Revt/10000
    
    if x < 10:
        x = 10
    
    index_vt   = nomes.index(perfil_EV)
    #Clalpha_vt = eval(clalphas[index_vt])
    Cd0_vt     = eval(cd0s[index_vt])
    #alpha0_vt  = eval(alphas_0[index_vt])
    t_c_vt     = t_cs[index_vt]
    
    math.sin(x) #so para não da erro chato, ignora

    
    '________________________________________________________________________________________'
    
    #asa
    
    b = lista_b_asa[-1]*2

    S = b * MAC_asa
    
    AR = b/MAC_asa
    
    CLmax,listax,listay = LLT(lista_b_asa,
                              lista_Cr_asa,
                              lista_Ct_asa,
                              Clalpha_w,
                              torção_asa,
                              a_s,
                              alpha0_w,
                              diedro_asa,
                              100)
    
    #EH
    
    b_EH = lista_b_EH[-1]*2

    S_EH = b_EH * MAC_EH
    
    AR_EH = b_EH / MAC_EH
        
    CLmax_EH,listax,listay = LLT(lista_b_EH,
                              lista_Cr_EH,
                              lista_Ct_EH,
                              Clalpha_EH,
                              0,
                              a_s_EH,
                              alpha0_EH,
                              0,
                              100)
    
    
    CLa_EH = (CLmax_EH - 0)/(a_s_EH - alpha0_EH)
    CL0_EH = CLmax_EH - (CLa_EH*a_s_EH)
    
    #Endplate

    M = V/343

    fM = 1 - 0.08*(M**1.45)

    Rend = V*Ced/vis

    CFed = 0.42/(np.log10(Rend)**2.58)

    ftced = 1 + (2.7*(t_end)) + (100*(t_end**4))

    CD0end = CFed * ftced * fM * (Swetend/S) * ((Cd0ed/0.004)**0.4)

    fator_area = -0.980392156862763*(Sep/S)**2 + 2.09803921568627*(Sep/S) + 0.999999999999998 # AR efetivo com endplate de acordo com a altura

    fator_altu = -0.980392156862763*(h_end/b)**2 + 2.09803921568627*(h_end/b) + 0.999999999999998 # AR efetivo com endplate de acordo com a altura

    fator = (fator_area+fator_altu)/2
    
    Ae = AR*fator
    
    e=1/(1.05+(0.007*np.pi*AR)) #oswald.pdf
 
    delt_CL = (np.pi*e*Ae*(2*CD0end*(Sep/S)))**0.5
    
    CLmax = CLmax + delt_CL #compensação do endplate

    CLa = (CLmax - 0)/(a_s - alpha0_w)
    CL0 = CLmax - (CLa*a_s)
        


    
    '________________________________________________________________________________________'
    
    M = V/343

    CFw  = 0.42/(np.log10(Re)**2.58)
    CFht = 0.42/(np.log10(Reht)**2.58)
    CFvt = 0.42/(np.log10(Revt)**2.58)
    CFf  = 0.42/(np.log10(Ref)**2.58)
    CFed = 0.42/(np.log10(Rend)**2.58)

    fM = 1 - 0.08*(M**1.45)

    fld = 1 + (60/((lf/df)**3))+(0.0025*(lf/df))

    ftcw  = 1 + (2.7*(t_c_w))  + (100*(t_c_w**4))
    ftcht = 1 + (2.7*(t_c_EH)) + (100*(t_c_EH**4))
    ftcvt = 1 + (2.7*(t_c_vt)) + (100*(t_c_vt**4))
    ftced = 1 + (2.7*(t_end)) + (100*(t_end**4))
    
    CD0w   = CFw  * ftcw  * fM * 2 * ((Cd0_w/0.004)**0.4)
    CD0ht  = CFht * ftcht * fM * 2 * ((Cd0_EH/0.004)**0.4)
    CD0vt  = CFvt * ftcvt * fM * 2 * ((Cd0_vt/0.004)**0.4)
    CD0f   = CFf  * fld   * fM * 4 * ((Cd0_f/0.004)**0.4)
    CD0end = CFed * ftced * fM * 2 * ((Cd0ed/0.004)**0.4)
    
    
    
    # SUSTENTACAO
    
    #ASA
    
    solo_CL_asa = efeito_solo(MAC_asa,b,h,V)
    
    CLa = CLa * solo_CL_asa #efeito solo
    
    CL = (CLa * (i_w + alpha)) + CL0
    
    L = (rho * (V**2) * S * CL)/2
    
    #EH
    
    solo_CL_EH = efeito_solo(MAC_EH,b_EH,h,V)
    
    CLa_EH = CLa_EH * solo_CL_EH #efeito solo
    
    CL_EH = (CLa_EH * (i_EH + alpha)) + CL0_EH
    
    L_EH = (rho * (V**2) * S_EH * CL_EH)/2
    

    L_total = L + L_EH
    
    
    # ARRASTO
    
    solo = ((16*h/b)**2)/(1+((16*h/b)**2))
    
    #ASA
    
    CD_ASA = solo*(CL**2)/(np.pi*e*Ae) + CD0w + 2*CD0end*(Sep/S)
    
    D = (rho * (V**2) * S * CD_ASA)/2
    
    #EH
    
    CD_EH = solo*(CL_EH**2)/(np.pi*e*AR_EH) + CD0ht
    
    D_EH = (rho * (V**2) * S_EH * CD_EH)/2
    
    #EV
    
    CD_EV =  CD0vt
    
    D_EV = (rho * (V**2) * S_EV * CD_EV)/2
    
    #Aproximação da fuselagem
    
    D_f = (rho * (V**2) * (lf*df) * CD0f)/2

    D_total = D + D_EH + D_EV + D_f
    
    #print()
    #print(V)
    #print(CL0,CLa,alpha,CL,)
    #print(L_total)
    # print(CD_ASA)
    
    return(L_total,D_total)



def remove_lista_da_outra(lista_original,lista_remover):
    
    for i in lista_remover:
        
        lista_original.remove(i)
    
    return(lista_original)

def efeito_solo(MAC,b,h,V):
    '_______________________________________________________________________________________'

    #Efeito do solo no CL (Nicolai pag 259)
    
    didiv_AR = + 0.22588727030313294 -1.7195700699150223*(2*h/b)*((2*h/b)-1) -0.0014789585540273542*(2*h/b)**(-1) + 0.8814936755669294*(2*h/b)**3 -0.181171052815559*(2*h/b)**4
    
    if 2*h/b > 1.8:
        
        didiv_AR = 1

    M = V/343

    beta = (1 - M**2)**0.5

    AR = b/MAC

    CLa = (2*np.pi*AR)/(2+(4+(AR**2 * beta**2 *(1+(0/beta**2))))**0.5)

    AR = AR/didiv_AR

    CLa2 = (2*np.pi*AR)/(2+(4+(AR**2 * beta**2 *(1+(0/beta**2))))**0.5)

    solo_CL_a = CLa2/CLa
    
    return(solo_CL_a)

    '_______________________________________________________________________________________'





xs = []
vs = []
ms = []
ts = []
subiu = []
ys = []
alphas = []
dx = []
dy = []

guarda_x = []
guarda_V = []

guarda_L = []
guarda_D = []

guarda_t = []
guarda_a = []

guarda_Vx = []
guarda_Vy = []

guarda_alpha = []


angulo_subida = []

print("Carregando...")

for m in np.arange(comecar_MTOW,MTOW_final,incremento_MTOW):
    
    print(f"\nTestando MTOW = {round(m,3)} Kg")
    
    Vx = 0 #velocidade x
    Vy = 0 #velocidade y
    V_t = 0 #resultante das velocidades x e y

    x = 0 #posição x
    y = 0 #posição y

    tempo = 0 # tempo 

    ax = 0  # aceleração em x
    ay = 0  # aceleração em y
    
    alpha = 0 #angulo do avião
    
    dina_x = [0] #começa no zero(guarda posições)
    dina_y = [0] #começa no zero(guarda posições)
    dina_t = [0] #começa no zero(guarda tempo)
    dina_a = [0] #começa no zero(guarda aceleração)
    dina_alpha = [0] #começa no zero(guarda alpha)
    dina_Vx = [0]  #começa no zero(guarda Vx)
    dina_Vy = [0]  #começa no zero(guarda Vy)
    
    g_x = 0 #valor só para guardar qunado ele saiu do chão
    g_v = 0 #valor só para guardar qunado ele saiu do chão
    
    L = 0
    D = 0
    
    memoria_L = [0]
    memoria_D = [0]
    
    
    while y < h_obs  and x <= pista_max:
        
        "O jeito que funciona é:"
        "Se a somatória de forças em y der positiva = avião decolou, para o loop"
        "Se o x passar do tamanho da pista, para o loop"
        
        "Avião no começo da pista, x = 0, V = 0 e a = 0"
        "0) usa velocidade do loop anterior(caso exista)"
        "1) calcula a somatoria de forças e em X"
        "2) calcula aceleração naquele momento, considerando a massa"
        "3) guarda a posicao de X inicial para calcular a velocidade dps"
        "4) calcula o novo X somando ao anterior e calculando com a equação de movimento linear-"
        "-contanto que o incremento de tempo seja pequeno o suficiente(<0.1) esse novo valor de X é confiável"
        "5) calcula o novo tempo"
        "6) calcula a nova velocidade"
        
        #A somatoria de forças em x é:
        #Tração do motor * cos(AOA) - atrito do chão(em função da força normal, que é PESO - L) - Arrasto aerodinâmico
        
        alpha_rad = alpha*np.pi/180
        
        tracao = eval(str(T))
        
        if Vx > 0:
            Fx = (np.cos(alpha_rad) * tracao) - atrito*((m*g)-(L*np.cos(alpha_rad))) - D*np.cos(alpha_rad)
        else:
            Fx = (np.cos(alpha_rad) * tracao) #se ele ta parado, só tem a tração
        #A somatoria de forças em y é:
        #Sustentação - Peso + Tração do motor * sen(AOA)
        
        Fy = L*np.cos(alpha_rad) - (g*m) + (np.sin(alpha_rad)*tracao) - (D*np.sin(alpha_rad))
        
        
        if Fy > 0: #se saiu do chão, não tem mais atrito das rodas
            Fx = (np.cos(alpha_rad) * tracao) - D*np.cos(alpha_rad)
        
        if Fy > 0 and g_x == 0 and g_v == 0: #isso aqui  é só para guardar quando ele decolou do chão
            g_x = x
            g_v = Vx
        
        
        ax = Fx/m #Eq de Newton
        ay = Fy/m #Eq de Newton
        
        dina_a.append(ax)
        
        x_anterior = x #guardando o valor de x antes de calcular o próximo
        y_anterior = y #guardando o valor de y antes de calcular o próximo
        
        x = x  + Vx*incremento_tempo + (ax*incremento_tempo**2)/2 #equação de movimento linear
        y = y  + Vy*incremento_tempo + (ay*incremento_tempo**2)/2 #equação de movimento linear
        
        
        if y < 0: #não tem como o avião está abaixo do chão (restrição de contorno)
            y = 0
            
        h = h_asa + y #a altura que será usada no eveito solo, é a altura da asa mais a posição em y
            
        dina_x.append(x) #guardando a posição para fazer animação no final
        dina_y.append(y) #guardando a posição para fazer animação no final

        tempo = tempo + incremento_tempo #proxima iteracao de tempo
        
        dina_t.append(tempo)
        
        Vx = (x - x_anterior)/incremento_tempo # nova velocidade, considerando quanto ele andou no tempo
        Vy = (y - y_anterior)/incremento_tempo # nova velocidade, considerando quanto ele subiu no tempo
        
        dina_Vx.append(Vx) #guardando a Vx para fazer grafico no final
        dina_Vy.append(Vy) #guardando a Vy para fazer grafico no final
        
        #a velocidade que será referente aos cálculos de L e D será a velocidade resultante
        V_t = (Vx**2 + Vy**2)**0.5 
        
        
        #O Alpha ou AOA ou ângulo de ataque é referente a resultante das velocidade dos avião, o ângulo da resultante

        #lembrar que arctan da em rad, tem q transformar em grau para somar com i_w

        
        if analise == "ideal" and x > x_ativa_profundor:         #caso seja ideal:
            #o profundor conseguiu rotacionar o avião para o CLmax da asa em 1 segundo
            alpha = alpha + (alpha_s - i_w)*(incremento_tempo) + np.arctan(Vy/Vx)*(180/np.pi)
        else:
            alpha = np.arctan(Vy/Vx)*(180/np.pi) 
            

        if alpha >= alpha_s - i_w: #quando o angulo de estol é alcançado, o angulo n sobe mais
            alpha = alpha_s - i_w
        
        
        dina_alpha.append(alpha)
         
        
        
        L, D = L_e_D(arquivos_perfil,
                    #asa
                    perfil_asa,
                    lista_b_asa,
                    lista_Cr_asa,
                    lista_Ct_asa,
                    torção_asa,
                    diedro_asa,
                    i_w,
                    h,
                    alpha,
                    
                    #EH
                    perfil_EH,
                    lista_b_EH,
                    lista_Cr_EH,
                    lista_Ct_EH,
                    
                    #EV
                    perfil_EV,
                    lista_b_EV,
                    lista_Cr_EV,
                    lista_Ct_EV,
                    
                    #Fuselagem
                    lf,
                    df,
                    
                    #Constantes
                    V_t,
                    rho,
                    vis,)
        
        L = L*np.cos(alpha*np.pi/180) - D*np.sin(alpha*np.pi/180)
        
        memoria_L.append(L)
        memoria_D.append(D)
        
        #aqui acaba o while

    
    "isso aqui é para ele parar caso os aviões estejam passando da pista sem subir pra voo"
    "Se a massa for muito grande ele n consegue dar voo no tamanho da pista, então a condição-"
    "-de x <= pista é alcançada sem a outra condição de Voar ser alcançada"
    if y < h_obs or x >= pista_max: 
        print("Não decolou!")
        break
        
    #tg(angulo_subida) = cateto oposto(y) / cateto adjacente(x)
    alpha_subida = np.arctan(h_obs/(pista_max - g_x)) * 180/np.pi
    
    print(f"Decolou em {round(g_x,2)} m e em {round(tempo,2)} segundos")
    print(f'Velocidade: {round(Vx,2)}')
    print(f"Ângulo de subida {round(alpha_subida)}°")
    print(f'AOA alcançado: {round(alpha,2)}°')
    
    
    
    angulo_subida.append(alpha_subida)
    
    ts.append(tempo)
    xs.append(x)
    vs.append(Vx)
    ms.append(m)
    ys.append(y)
    alphas.append(alpha)
    #subiu.append(eval(somatorio_Fy))
    
    dx.append(dina_x)
    dy.append(dina_y)
    guarda_t.append(dina_t)
    guarda_a.append(dina_a)
    
    guarda_x.append(g_x)
    guarda_V.append(g_v)
    
    guarda_L.append(memoria_L)
    guarda_D.append(memoria_D)
    
    guarda_Vx.append(dina_Vx)
    guarda_Vy.append(dina_Vy)
    
    guarda_alpha.append(dina_alpha)


"aqui eu pego o ultimo valor das listas que eu guardei pq vai ser os valores máximos, ou seja do MTOW"

MTOW = round(ms[-1],2)                
V_decolagem = round(vs[-1],2)
x_pista = round(xs[-1],2)
t_decolagem = round(ts[-1],2)
altura_aviao_no_x_pista = round(ys[-1],2)
apha_alcancado = round(alphas[-1],2)
angulo_de_subida = round(angulo_subida[-1],2)

x_que_subiu = round(guarda_x[-1],2)
V_que_subiu = round(guarda_V[-1],2)

L_durante_decolagem = guarda_L[-1]
D_durante_decolagem = guarda_D[-1]

Tempo_decolagem_dinamico = guarda_t[-1]
Aceleracao_decolagem_dinamico = guarda_a[-1]

Vx_decolagem_dinamico = guarda_Vx[-1]
Vy_decolagem_dinamico = guarda_Vy[-1]

vx = dx[-1]
vy = dy[-1]

alpha_dinamico = guarda_alpha[-1]


aceleracao_max = round(max(Aceleracao_decolagem_dinamico),2)

RC = round(integral_sem_integral(Tempo_decolagem_dinamico, vy),2)


print("___________________________________________")
print(f"\n{MTOW = } Kg ({analise})")
print(f"{V_decolagem = } m/s")
print(f"{x_pista = } m")
print(f"{altura_aviao_no_x_pista = } m")
print(f"{t_decolagem = } s")
print(f"{aceleracao_max = } m/s²")
print(f"{angulo_de_subida = }°")
print(f"{RC = } m/s")

print(f"\n{x_que_subiu = } m")
print(f"{V_que_subiu = } m/s")
print(f"{apha_alcancado = }°")
print("___________________________________________")

try:
    plt.style.use('extensys-gd')
except:
    pass



plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(vx,Vx_decolagem_dinamico,color = "blue")
plt.axvline(x_que_subiu,label="Decolagem",color = "black")
plt.title("Vx vs x")
plt.xlabel("[m] x")
plt.ylabel("[m/s] Vx")
plt.legend(bbox_to_anchor=(1.02, 0.85, 0.25, 0.1), frameon=True, shadow=True, ncol=1)
plt.show()

plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(ms,ts,color = "blue")
plt.title("Tempo da decolagem em função da massa do avião")
plt.xlabel("[Kg] Massa total do avião")
plt.ylabel("[s] Tempo para decolagem")
plt.show()


plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(ms,xs,color = "green")
plt.title("Pista de decolagem necessária em função da massa do avião")
plt.xlabel("[Kg] Massa total do avião")
plt.ylabel("[m] Pista necessária para decolagem")
plt.show()


plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(vx,vy,color = "blue",label = "Avião")
plt.title("Trajetória do avião")
plt.xlabel("[m] Posição x do avião")
plt.ylabel("[m] Altura")
plt.axis("scaled")
plt.axhline(-0.1,color= "grey",label = "Chão")
plt.plot([pista_max,pista_max],[0,h_obs],color = "red",label = "Obstaculo")
plt.axis([0,60,-1,3])
plt.legend(bbox_to_anchor=(1.4,0.6), frameon=True, shadow=True, ncol=2)
plt.show()

#graficos de aerodinamica

# L_D = []

# for i in range(len(vx)):
    
#     if D_durante_decolagem[i] == 0:
#         L_D.append(0)
#     else:
#         L_D.append(L_durante_decolagem[i]/D_durante_decolagem[i])

# plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
# plt.plot(vx,L_durante_decolagem,color = "blue",label="L")
# plt.title("Sustentação gerada pelo avião em Posição na pista")
# plt.ylabel("[N] Sustentação gerada pelo avião")
# plt.xlabel("[m] Posição na pista")
# plt.axvline(x_que_subiu,label="Decolagem",color = "black")
# plt.legend(bbox_to_anchor=(1.02, 0.85, 0.25, 0.1), frameon=True, shadow=True, ncol=1)
# plt.show()

# plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
# plt.plot(vx,D_durante_decolagem,color = "red",label="D")
# plt.title("Arrasto gerada pelo avião em Posição na pista")
# plt.ylabel("[N] Arrasto gerado pelo avião")
# plt.xlabel("[m] Posição na pista")
# plt.axvline(x_que_subiu,label="Decolagem",color = "black")
# plt.legend(bbox_to_anchor=(1.02, 0.85, 0.25, 0.1), frameon=True, shadow=True, ncol=1)
# plt.show()

# plt.figure(num=None, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
# plt.plot(vx,L_D,color = "green",label="L/D")
# plt.title("Eficiência do avião em Posição na pista")
# plt.ylabel("L/D")
# plt.xlabel("[m] Posição na pista")
# plt.axvline(x_que_subiu,label="Decolagem",color = "black")
# plt.legend(bbox_to_anchor=(1.02, 0.85, 0.25, 0.1), frameon=True, shadow=True, ncol=1)
# plt.show()


#animação

guarda_ultimo_ponto_x = vx[-1]
guarda_ultimo_ponto_y = vy[-1]

for i in range(4):
    
    removerx = []
    removery = []
    
    for i in range(len(vx)):
        
        if i % 2 == 1:
            removerx.append(vx[i])
            removery.append(vy[i])
            
    
    
    vx = remove_lista_da_outra(vx,removerx)
    vy = remove_lista_da_outra(vy,removery)

vx.append(guarda_ultimo_ponto_x)
vy.append(guarda_ultimo_ponto_y)


fig, ax = plt.subplots()
scat = ax.scatter(vx[0], vy[0], c="b",s = 5)
ax.axhline(-0.1,color= "grey",label = "Chão")
ax.plot([pista_max,pista_max],[0,h_obs],color = "red",label = "Obstaculo")
ax.axis("scaled")
ax.axis([0,60,-1,3])
#plt.legend(bbox_to_anchor=(1.3,0.15), frameon=True, shadow=True, ncol=2)



def update(frame):
    # for each frame, update the data stored on each artist.
    x = vx[:frame]
    y = vy[:frame]
    # update the scatter plot:
    data = np.stack([x, y]).T
    scat.set_offsets(data)
    return (scat)

print("\nFazendo animação...")

mili_seg_por_ponto = round(1000 * t_decolagem/len(vx))

#interval = 1000 #significa 1 frame a cada 1 segundo

ani = animation.FuncAnimation(fig=fig, func=update, frames = len(vx), interval=int(mili_seg_por_ponto))


ani.save(filename="Tragetoria aviao.gif", writer="pillow")

print("Animação feita!")









