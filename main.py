import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import geometria
import SistemaRD_Schnakenberg
import homogeneizacion
from scipy.optimize._numdiff import approx_derivative
import math

#No Borrar esta parte, es para guardar el codigo en la carpeta de resultados
import os
import shutil
# crear carpeta si no existe
os.makedirs("Resultados", exist_ok=True)
# nombre del archivo actual
archivo_origen = "main.py"
archivo_destino = f"Resultados/main.py"
# copiar archivo
shutil.copy(archivo_origen, archivo_destino)
print("Archivo guardado en:", archivo_destino)
carpera_resultados="Resultados"


# -------------------------------------------------
# Definiendo la geomertría del modelo
# -------------------------------------------------
# Parametros
Sfin=2*np.pi
r=0.3
R=1
eps=1/(7)
p=40
q=5/eps   
#Construct
p_geometria=[R,r,p,q,Sfin]
geomToro=geometria.geometria_class(1,p_geometria)
Sfin=geomToro.p_geom[-1]
print("Sfin",Sfin)

Grafica_Prop_Geometricas=True
if Grafica_Prop_Geometricas:
    geomToro.grafica_curva(0,Sfin)
    geomToro.GraficaCurvaura(0,Sfin,FileName="Resultados/Curvatura")
    #geomToro.GraficaVTangente(0,Sfin)
    # print(geomToro.LongitudArco(0,np.pi/4))


#Esfera
# R=1.0
# eps=1/10.0
# a=1.0/eps
# b=1.0/eps
# p_geometria=[R,a,b]
# geomEsf=geometria.geometria_class(2,p_geometria)
# geomEsf.grafica_curva(0,Sfin)
# geomEsf.imprimir_mi_tipo()

# ------------------------------------------------------
# Deficnicion de el sistema sobre la curca
# ------------------------------------------------------
# Parametros
# 
a=0.1
b=1.5
d=0.052
eta=0.15
alpha=0.01   #modelo (1) de difusión con curvatura seccion (4.1)
p=[a,b,d,eta,alpha] 

SisRD=SistemaRD_Schnakenberg.RD_Class(p,geomToro)
SisRD.GraficaDifusion(0,Sfin,num_points=10000)


#
# ------------------------------------------------------
# Definicion Homogeneización
# ------------------------------------------------------
# Parametros
# 
H_RD=homogeneizacion.H_RD_Class(SisRD)

promu,promv=H_RD.Coef_Efectivo(0,2*np.pi)
print("Coeficiente de difusión efectivo para u:", promu)
print("Coeficiente de difusión efectivo para v:", promv)

#Analisis de inestabilidad de Turing medio heterogéneo
InestabilidadTuringHet=False
if InestabilidadTuringHet:
    d_min=0.02
    d_max=0.06
    alpha_min=0.05
    alpha_max=1.5
    N1=100
    inter_d=np.linspace(d_min, d_max, N1)
    inter_alpha=np.linspace(alpha_min, alpha_max, N1)
    H_RD.InestabilidadHomogeneizado(inter_d,inter_alpha,labelx="d",labely="alpha")


# -------------------------------------------------
# Simulaciones
# -------------------------------------------------
Simulacion= True
times_to_save0=[0,10,20,30,40,50,100,150,200,500,600,700,800,900,1000]
if Simulacion:
    print("------------Simulación con tiempos de salida---------------------")
    print("times_to_save:",times_to_save0)
    
    H_RD.SimulHomogenizacion(p,S=SisRD.geom.p_geom[-1], num_points=600, times_to_save=times_to_save0, periodic=True,FileName="Resultados/solucionHomo",dte=1e-5)
    H_RD.SimulRD(p,S=SisRD.geom.p_geom[-1], num_points=600, times_to_save=times_to_save0, periodic=True,FileName="Resultados/solucionaHetero",dte=1e-5 )

# -------------------------------------------------
# Grafica Comparacion de Soluciones
# -------------------------------------------------
Grafica_CompararSol=True
if Grafica_CompararSol:
    NombreDelArchivoU=f"{carpera_resultados}/solucionaHetero_1000.00.txt"
    NombreDelArchivoH=f"{carpera_resultados}/solucionHomo_1000.00.txt"
    Salida=f"{carpera_resultados}/Comparacion1000.png"
    H_RD.GrafComparacion(FilenameU=NombreDelArchivoU,FilenameH=NombreDelArchivoH,FileSalida=Salida)

# -------------------------------------------------
# Animacion de Soluciones
# -------------------------------------------------
Animacion_Sol=False
if Animacion_Sol:
    geomToro.generar_gif(1, "u(s,t)", filename="Resultados/solucionaHetero", salida_gif="Resultados/animacionHetero.gif",tiempos=times_to_save0)

# -------------------------------------------------
# Grafica 3D de Archivo Guradado en txt: #s f
# -------------------------------------------------
Grfica_ArchivoTxt=True
if Grfica_ArchivoTxt:
    NombreDelArchivo=f"{carpera_resultados}/solucionaHetero_1000.00.txt"
    NombreSalida=f"{carpera_resultados}/solucionaHetero_1000.00.png"
    NombreFuncion0="u(s,t=1000)"
    geomToro.GrafColorF(Filename=NombreDelArchivo,FileSalida=NombreSalida,NombreFuncion=NombreFuncion0,barra=False)

# -------------------------------------------------
# Grafica 3D Curvatura de Archivo Guradado en txt: #s f
# -------------------------------------------------
Grfica_ArchivoTxt=True
if Grfica_ArchivoTxt:
    NombreDelArchivo=f"{carpera_resultados}/Curvatura.txt"
    NombreSalida=f"{carpera_resultados}/Curvatura3D.png"
    NombreFuncion0="kg(s)"
    geomToro.GrafColorF(Filename=NombreDelArchivo,FileSalida=NombreSalida,NombreFuncion=NombreFuncion0)


Grafica_ArchivoTxtProyeccion=False
if Grafica_ArchivoTxtProyeccion:
    NombreDelArchivo=f"{carpera_resultados}/solucionaHetero_1000.00.txt"
    NombreSalida=f"{carpera_resultados}/ProySolucionaHetero_1000.00.png"
    NombreFuncion0="u(s,t=1000)"
    geomToro.GrafColorF(Filename=NombreDelArchivo,FileSalida=NombreSalida,NombreFuncion=NombreFuncion0)
    geomToro.proyeccionSobreElPlano(Filename=NombreDelArchivo,FileSalida=NombreSalida,NombreFuncion=NombreFuncion0)
    #geomToro.proyeccionSobreElPlanoScatter(Filename=NombreDelArchivo,FileSalida=NombreSalida,NombreFuncion=NombreFuncion0)