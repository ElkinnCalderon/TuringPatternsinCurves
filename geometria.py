import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize._numdiff import approx_derivative
import math
import glob
import imageio
import os

class geometria_class:
    def __init__(self, tipo,p_geom):
        #self.L=L
        self.curva=tipo  #1: toro, 2:esfera,....
        self.p_geom=p_geom  #parámetros de la curva, por ejemplo (R,r,p,q) para el toro
        self.cerrada=False

        if self.curva==1:
            from fractions import Fraction
            from math import gcd

            p = self.p_geom[2]
            q = self.p_geom[3]

            print("vueltaToroidal p:", 2*np.pi / p)
            print("vueltaPoloidal q:", 2*np.pi / q)

            razon = p / q
            print("p/q =", razon)

            # Aproximar p/q como fracción racional
            frac = Fraction(razon).limit_denominator(100)

            m = frac.numerator
            n = frac.denominator

            print(f"p/q ≈ {m}/{n}")

            # Si p/q es racional, la curva se cierra
            if np.isclose(razon, m / n):
                print("La curva se cierra y vuelve al punto inicial.")
                self.cerrada=True

                # valor mínimo de s tal que:
                # p*s = 2π M
                # q*s = 2π N

                # usando:
                # s = 2π * n / q = 2π * m / p

                s_cierre = 2 * np.pi * n / q

                print("Valor mínimo de s donde vuelve al inicio:", s_cierre)
                if self.p_geom[-1] != s_cierre:
                    
                    self.p_geom[-1] = s_cierre
                    print("Se modifico Sfinal, Nuevo:",self.p_geom[-1])
                    

            else:
                print("La curva no vuelve exactamente.")
                print("Rellena densamente el toro.")

    def curvatura(self,t):
        if self.curva==1:
            R = self.p_geom[0]  # Radio mayor del toro
            r = self.p_geom[1]  # Radio menor del toro
            p=self.p_geom[2]
            q=self.p_geom[3]
            u = p * t
            v = q * t

            # Curvatura de una curva sobre el toro
            kn= -np.cos(v) / (p+ q * np.cos(v))

            #C
            numerador = p * np.sin(q*t) * (p**2 * (R + r*np.cos(q*t))**2 + 2 * r**2 * q**2)
            denominador = (p**2 * (R + r*np.cos(q*t))**2 + r**2 * q**2)**(3/2)
            kg = numerador / denominador
            
            k=kn+kg
        elif self.curva==2:
            R = self.p_geom[0]  # Radio de la esfera
            a=self.p_geom[1]
            eps=self.p_geom[2]

            # Curvatura de una curva sobre la esfera
            k = (a**2 + 1) / (R * np.sin(a**(t/eps)))
            kg = np.sqrt(((a**2 + 1)**2) / (R**2 * np.sin(a**(t/eps))**2) - 1/(R**2))
            kg= np.sqrt(k**2-kg**2)

        return k,kg,kn
    
    def GraficaCurvaura(self, a, b,ik=1, num_points=1000,FileName="Curvatura"):
        ffig, ax1 = plt.subplots(
            figsize=(8, 6)
        )

        intT = np.linspace(a, b, num_points)
        kCurvatura = []
        for t in intT:
            ks= self.curvatura(t)
            kCurvatura.append(ks[ik])
        
        nombre_archivo = f"{FileName}.txt"
        datos = np.column_stack((intT, kCurvatura))

        np.savetxt(
                nombre_archivo,
                datos,
                header="x u v",
                fmt="%.8f"
            )

        ax1.plot(intT, kCurvatura, label='Curvatura')
        #ax1.axhline(dc*0.0002773314352045536, linestyle="--")
        ax1.set_ylabel("k")
        ax1.set_title("Curvatura")
        ax1.grid(True)


        plt.tight_layout()
        plt.show()
    
    def LongitudArco(self, a, b, num_points=1000):
        intT = np.linspace(a, b, num_points)
        longitud = 0.0
        longitud2= 0.0

        
        for i in range(1, len(intT)):
            t1 = intT[i-1]
            t2 = intT[i]
            DeltaT=t2-t1
            x1,y1,z1,ds = self.VectorTangente((t1+t2)/2)
            ds=ds*DeltaT
            # Aproximación de la longitud de arco usando la fórmula de la longitud de arco derivada
            longitud += ds

            x1,y1,z1 = self.FvectorialCurva(t1)
            x2,y2,z2 = self.FvectorialCurva(t2)
            # Aproximación de la longitud de arco usando la fórmula de la longitud de arco derivada
            ds2=np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            longitud2 += ds2


        return longitud,longitud2
    

    def FvectorialCurva(self,s):
        # Ángulos
        if self.curva==1:
            R = self.p_geom[0]  # Radio mayor del toro
            r = self.p_geom[1]  # Radio menor del toro
            u = self.p_geom[2] * s
            v = self.p_geom[3] * s

            # Parametrización del toro
            x = (R + r * np.cos(v)) * np.cos(u)  #r * -self.p_geom[3] np.sin(v)np.cos(u)-(R + r * np.cos(v)) * self.p_geom[2] np.sin(u)
            y = (R + r * np.cos(v)) * np.sin(u) #r * -self.p_geom[3] np.sin(v)np.sin(u)+(R + r * np.cos(v)) * self.p_geom[2] np.cos(u)
            z = r * np.sin(v) #r * np.cos(v) * self.p_geom[3]

        elif self.curva==2:
            # Esfera  # Radio de la esfera
            R = self.p_geom[0]  # Radio de la esfera
            a=self.p_geom[1]
            b=self.p_geom[2]
            u =  a*s
            v =  b*s

            # Parametrización de la esfera
            x = R * np.sin(v) * np.cos(u) #R * (a/eps) * np.cos(v) * np.cos(u) - R * (1/eps) * np.sin(v) * np.sin(u)
            y = R * np.sin(v) * np.sin(u) #R * (a/eps) * np.cos(v) * np.sin(u) + R * (1/eps) * np.sin(v) * np.cos(u)
            z = R * np.cos(v)#- R * (a/eps) * np.sin(v)
        
        return x, y, z
    
    def VectorTangente(self,s):
        # Cálculo del vector tangente a la curva
        if self.curva==1:
            R = self.p_geom[0]  # Radio mayor del toro
            r = self.p_geom[1]  # Radio menor del toro
            u = self.p_geom[2] * s
            v = self.p_geom[3] * s

            # Derivada de la parametrización del toro
            dx_dt = r * -self.p_geom[3] * np.sin(v) * np.cos(u) - (R + r * np.cos(v)) * self.p_geom[2] * np.sin(u)
            dy_dt = r * -self.p_geom[3] * np.sin(v) * np.sin(u) + (R + r * np.cos(v)) * self.p_geom[2] * np.cos(u)
            dz_dt = r * np.cos(v) * self.p_geom[3]

        elif self.curva==2:
            R = self.p_geom[0]  # Radio de la esfera
            a=self.p_geom[1]
            eps=self.p_geom[2]
            u =  s
            v =  a * s/eps

            # Derivada de la parametrización de la esfera
            dx_dt = R * (a/eps) * np.cos(v) * np.cos(u) - R * (1/eps) * np.sin(v) * np.sin(u)
            dy_dt = R * (a/eps) * np.cos(v) * np.sin(u) + R * (1/eps) * np.sin(v) * np.cos(u)
            dz_dt = -R * (a/eps) * np.sin(v)

        return dx_dt, dy_dt, dz_dt, np.sqrt((dx_dt)**2 + (dy_dt)**2 + (dz_dt)**2)
    
    def GraficaVTangente(self, a, b,ik=3, num_points=1000):
        ffig, ax1 = plt.subplots(
            figsize=(8, 6)
        )

        intT = np.linspace(a, b, num_points)
        kCurvatura = []
        for t in intT:
            ks= self.VectorTangente(t)
            kCurvatura.append(ks[ik])

        ax1.plot(intT, kCurvatura, label='Difusión de u')
        #ax1.axhline(dc*0.0002773314352045536, linestyle="--")
        ax1.set_ylabel("k")
        ax1.set_title("Vector Tangente")
        ax1.grid(True)


        plt.tight_layout()
        plt.show()
    
    def grafica_curva(self, a, b, num_points=10000):
        intT = np.linspace(a, b, num_points)

        x, y, z = [], [], []
        for t in intT:
            xi, yi, zi = self.FvectorialCurva(t)
            x.append(xi)
            y.append(yi)
            z.append(zi)

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        # --- Figura grande ---
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # --- Gráfica ---
        ax.plot(x, y, z, linewidth=1.5)

        # --- Escala uniforme (clave) ---
        max_range = np.array([
            x.max() - x.min(),
            y.max() - y.min(),
            z.max() - z.min()
        ]).max() / 2.0

        mid_x = (x.max() + x.min()) / 2.0
        mid_y = (y.max() + y.min()) / 2.0
        mid_z = (z.max() + z.min()) / 2.0

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        ax.set_box_aspect([1, 1, 1])  # proporción real

        # --- Etiquetas ---
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Curva en el espacio')

        # --- Expandir el gráfico dentro de la figura ---
        ax.set_position([0.05, 0.05, 0.9, 0.9])

        # --- Opcional: quitar rejilla para estilo limpio ---
        ax.grid(False)

        # --- Ajuste automático ---
        plt.tight_layout()

        plt.show()

    def promedioSobreCurva(self,a,b,num_points=1000):
        intT = np.linspace(a, b, num_points)
        k_values = []
        for t in intT:
            k,kg,kn = self.curvatura(t)
            k_values.append(k)

        promedio_k = np.mean(k_values)
        return promedio_k
    
    def imprimir_mi_tipo(self):
        if self.curva == 1:
            print("La curva es un toro.")
        elif self.curva == 2:
            print("La curva es una esfera.")
        else:
            print("Tipo de curva desconocido.")

    

    def GrafColorF(self,Filename=None,FileSalida="GrafColor.png",NombreFuncion="f(s)",ik=1,barra=True):
        
        if Filename is not None:
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            data = np.loadtxt(Filename)

            # columnas
            t = data[:, 0]
            f = data[:, ik]


            # -------------------------------------------------
            # construir segmentos 3D
            # -------------------------------------------------
            x,y,z=self.FvectorialCurva(t)
            points = np.array([x, y, z]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # -------------------------------------------------
            # colección de líneas coloreadas
            # -------------------------------------------------

            lc = Line3DCollection(
                segments,
                cmap="viridis",
                norm=plt.Normalize(f.min(), f.max())
            )

            lc.set_array(f[:-1])   # tamaño compatible
            lc.set_linewidth(2)

            # -------------------------------------------------
            # graficar
            # -------------------------------------------------
            
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection="3d")

            plt.tight_layout()

            ax.add_collection3d(lc)

            # ajustar límites manualmente
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(y.min(), y.max())
            ax.set_zlim(z.min(), z.max())

            max_range = np.array([
                x.max() - x.min(),
                y.max() - y.min(),
                z.max() - z.min()
            ]).max() / 2.0

            mid_x = (x.max() + x.min()) / 2.0
            mid_y = (y.max() + y.min()) / 2.0
            mid_z = (z.max() + z.min()) / 2.0

            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            ax.set_box_aspect([1, 1, 1])

            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_title(NombreFuncion)

            if barra:
                plt.colorbar(
                    lc,
                    ax=ax,
                    label=NombreFuncion,
                    shrink=0.5,     # altura
                    aspect=15,      # más grande = más delgada
                    pad=0.05        # separación del gráfico
                )

            # quitar rejilla
            ax.grid(False)

            # quitar ejes completos
            ax.set_axis_off()

            plt.savefig(FileSalida, dpi=300)

            plt.show()
            
        else:
            print("Tiene que crear una achivo .txt: #x f")
            
    def generar_gif(self, ik, NombreFuncion,filename="solucionaHetero",salida_gif="animacion.gif",tiempos=None):
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        import glob, os, imageio

        # -------------------------------------------------
        # 1. obtener archivos
        # -------------------------------------------------
        if tiempos is not None:
            archivos = []
            for t in tiempos:
                nombre_archivo = f"{filename}_{t:.2f}.txt"
                if os.path.exists(nombre_archivo):
                    archivos.append(nombre_archivo)
                else:
                    print(f"Archivo no encontrado: {nombre_archivo}")
        else:
            archivos = sorted(glob.glob(f"{filename}*.txt"))

        if len(archivos) == 0:
            print("No se encontraron archivos")
            return

        imagenes = []

        # -------------------------------------------------
        # normalización global de color
        # -------------------------------------------------
        f_global_min = +np.inf
        f_global_max = -np.inf

        for Filename in archivos:
            data = np.loadtxt(Filename)
            f = data[:, ik]
            f_global_min = min(f_global_min, f.min())
            f_global_max = max(f_global_max, f.max())

        norm_global = plt.Normalize(f_global_min, f_global_max)

        # -------------------------------------------------
        # 2. loop principal
        # -------------------------------------------------
        for i, Filename in enumerate(archivos):

            data = np.loadtxt(Filename)

            t = data[:, 0]
            f = data[:, ik]

            x, y, z = self.FvectorialCurva(t)

            # segmentos
            points = np.array([x, y, z]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            #DOS colecciones (clave)
            lc1 = Line3DCollection(segments, cmap="viridis", norm=norm_global)
            lc2 = Line3DCollection(segments, cmap="viridis", norm=norm_global)

            lc1.set_array(f[:-1])
            lc2.set_array(f[:-1])

            lc1.set_linewidth(2)
            lc2.set_linewidth(2)

            # -------------------------------------------------
            # FIGURA con dos vistas
            # -------------------------------------------------
            fig = plt.figure(figsize=(16, 8))

            ax1 = fig.add_subplot(121, projection="3d")
            ax2 = fig.add_subplot(122, projection="3d")

            ax1.add_collection3d(lc1)
            ax2.add_collection3d(lc2)

            # vista superior
            ax2.view_init(elev=90, azim=0)

            # -------------------------------------------------
            # límites iguales
            # -------------------------------------------------
            max_range = np.array([
                x.max() - x.min(),
                y.max() - y.min(),
                z.max() - z.min()
            ]).max() / 2.0

            mid_x = (x.max() + x.min()) / 2.0
            mid_y = (y.max() + y.min()) / 2.0
            mid_z = (z.max() + z.min()) / 2.0

            for ax in [ax1, ax2]:
                ax.set_xlim(mid_x - max_range, mid_x + max_range)
                ax.set_ylim(mid_y - max_range, mid_y + max_range)
                ax.set_zlim(mid_z - max_range, mid_z + max_range)
                ax.set_box_aspect([1, 1, 1])
                ax.set_axis_off()

            # títulos
            ax1.set_title(Filename)
            ax2.set_title("Vista superior")

            # -------------------------------------------------
            # colorbar SOLO en la derecha
            # -------------------------------------------------
            plt.colorbar(
                lc2,
                ax=ax2,
                label=NombreFuncion,
                shrink=0.6,
                aspect=20,
                pad=0.05
            )

            # -------------------------------------------------
            # guardar imagen
            # -------------------------------------------------
            base = os.path.splitext(os.path.basename(Filename))[0]
            nombre_img = f"Resultados/{base}.png"

            plt.savefig(nombre_img, dpi=150)
            plt.close(fig)

            imagenes.append(imageio.imread(nombre_img))

        # -------------------------------------------------
        # 3. crear GIF
        # -------------------------------------------------
        imageio.mimsave(salida_gif, imagenes, fps=5)

        print("GIF generado en:", salida_gif)

    def proyeccion(self, s):
        p = self.p_geom[2]
        q = self.p_geom[3]

        s = np.asarray(s)

        u = (p * s) % (2 * np.pi)
        v = (q * s) % (2 * np.pi)

        return u, v


    # -------------------------------------------------
    # gráfico en el plano (sin saltos)
    # -------------------------------------------------
    def proyeccionSobreElPlano(self, Filename=None, FileSalida="GrafColor.png", NombreFuncion="f(s)", ik=1):

        if Filename is None:
            print("Tiene que crear un archivo .txt: # t f")
            return
        from matplotlib.collections import LineCollection
        data = np.loadtxt(Filename)

        # columnas
        t = data[:, 0]
        f = data[:, ik]

        # -------------------------------------------------
        # proyección (u,v)
        # -------------------------------------------------
        u, v = self.proyeccion(t)

        # -------------------------------------------------
        # detectar saltos por periodicidad
        # -------------------------------------------------
        du = np.abs(np.diff(u))
        dv = np.abs(np.diff(v))

        umbral = np.pi
        mask = (du < umbral) & (dv < umbral)

        # -------------------------------------------------
        # construir segmentos válidos (sin saltos)
        # -------------------------------------------------
        points = np.array([u, v]).T

        p1 = points[:-1][mask]
        p2 = points[1:][mask]

        segments = np.stack([p1, p2], axis=1)
        colors = f[:-1][mask]

        # -------------------------------------------------
        # colección de líneas coloreadas
        # -------------------------------------------------
        lc = LineCollection(
            segments,
            cmap="viridis",
            norm=plt.Normalize(f.min(), f.max())
        )

        lc.set_array(colors)
        lc.set_linewidth(2)

        # -------------------------------------------------
        # graficar
        # -------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))

        ax.add_collection(lc)

        ax.set_xlim(0, 2*np.pi)
        ax.set_ylim(0, 2*np.pi)

        ax.set_xlabel("u")
        ax.set_ylabel("v")
        ax.set_title(NombreFuncion)

        ax.set_aspect('equal')

        plt.colorbar(lc, ax=ax, label=NombreFuncion)
        plt.savefig(FileSalida, dpi=300)
        plt.show()

    def proyeccionSobreElPlanoScatter(self, Filename=None, FileSalida="GrafColor.png", NombreFuncion="f(s)", ik=1):

        if Filename is None:
            print("Tiene que crear un archivo .txt: # t f")
            return

        data = np.loadtxt(Filename)

        # columnas
        t = data[:, 0]
        f = data[:, ik]

        # -------------------------------------------------
        # proyección (u,v)
        # -------------------------------------------------
        u, v = self.proyeccion(t)

        # -------------------------------------------------
        # graficar puntos coloreados
        # -------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))

        sc = ax.scatter(
            u, v,
            c=f,               # color según f
            cmap="viridis",
            s=20               # tamaño de los puntos
        )

        ax.set_xlim(u.min(), u.max())
        ax.set_ylim(v.min(), v.max())

        ax.set_xlabel("u")
        ax.set_ylabel("v")
        ax.set_title(NombreFuncion)

        ax.set_aspect('equal')

        plt.colorbar(sc, ax=ax, label=NombreFuncion)
        plt.savefig(FileSalida, dpi=300)
        plt.show()