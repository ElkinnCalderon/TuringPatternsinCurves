import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import geometria
from scipy.optimize._numdiff import approx_derivative
import math

class RD_Class:
    def __init__(self, p,geom):
        self.p= p
        self.geom=geom
        print('Se definio una estructura sistema reaccion Difusion de Schnakenberg')

    def Reaccion(self, u, v):
        eta = self.p[3]
        a = self.p[0]
        b = self.p[1]
        f = eta * (a - u + (u**2.0) * v)
        g = eta * (b - (u**2.0) * v)
        return f, g

    def difusion(self, s,param1=None,param2=None):  ## Modelos definidos en la seccion 4.1 Tesis Dayana
        if param1==None and param2==None:
            d = self.p[2]
            alf=self.p[4]
        else:
            d=param1
            alf=param2
        k,kg,kn= self.geom.curvatura(s)
         #Modelo 1
        #du = d*(1+alf*kg**2)
        #dv = 1.0
        # Modelo 2
        du = d/(1+alf*k**2)
        dv = 1.0  
        #  Modelo 3, el parametro alpha en este  caso sera el umbral de curvatura dle modelo 3
        # d1=self.p[5]  # difusion para k>alf
        # if k<alf:
        #     du = d    
        #     dv = 1.0
        # else:     
        #     du = d1  ##### Importante: falta definirlo en main.py
        #     dv = 1.0
        return du, dv
    
    def GraficaDifusion(self, a, b, num_points=10000):
        fig, (ax1, ax2) = plt.subplots(
                2, 1,
                figsize=(8, 6),
                sharex=True
            )

        dc=self.dcritico()
        intT = np.linspace(a, b, num_points)
        difu_u = []
        difu_v = []
        for t in intT:
            du, dv = self.difusion(t)
            difu_u.append(du)
            difu_v.append(dv)

        ax1.plot(intT, difu_u, label='Difusión de u')
        #ax1.axhline(dc*0.0002773314352045536, linestyle="--")
        ax1.set_ylabel("Du")
        ax1.set_title(f"Coeficientes de Difusión")
        ax1.grid(True)


        ax2.plot(intT, difu_v, label='Difusión de v')
        #ax2.axhline(dc/0.00554669288555708, linestyle="--")
        ax2.set_xlabel("x")
        ax2.set_ylabel("Dv")
        ax2.grid(True)

        plt.tight_layout()
        plt.show()


    def Jacobiano(self,p, u, v):
        x0 = np.array([u, v], dtype=float)
        def F(x):
            f, g = self.Reaccion(x[0], x[1])
            return np.array([f, g])
        return approx_derivative(F, x0, method='2-point')

    def Equilibrio(self,p):
        """
        Calcula el punto de equilibrio (u_e, v_e)
        del sistema de Schnakenberg definido en Reaccion.
        """
        a = p[0]
        b = p[1]
        u_e = a + b
        v_e = b / (a + b)**2
        return u_e, v_e

    def Jacobiano_en_equilibrio(self,p):
        """
        Calcula el Jacobiano del sistema en el equilibrio (u_e, v_e)
        usando diferencias finitas.
        """
        u_e, v_e = self.Equilibrio(p)
        J = self.Jacobiano(p, u_e, v_e)
        return J, (u_e, v_e)

    def IntgracionEulerExplicito(p, u0, v0, dt, N):
        u = u0
        v = v0
        U = [u]
        V = [v]
        T = [0.0]
        for n in range(N):
            f, g = Reaccion(p, u, v)
            u = u + dt * f
            v = v + dt * g
            U.append(u)
            V.append(v)
            T.append((n+1)*dt)
        return np.array(U), np.array(V), np.array(T)
    def grafUVFase(U, V):
        plt.plot(U, V, 'o-', markersize=0.5, linewidth=1, label="Puntos")  # 'o-' dibuja puntos unidos por líneas
        plt.xlabel("u")
        plt.ylabel("v")
        plt.title("ListPlot en Python")
        plt.legend()
        plt.grid(True)
        plt.show()
    def grafUVTresD(U, V, T):
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(U, V, T, label='Trayectoria en 3D')
        ax.set_xlabel('u')
        ax.set_ylabel('v')
        ax.set_zlabel('t')
        ax.set_title('Gráfico 3D de u, v y t')
        plt.legend()
        plt.show()
    def grafUVTemporal(U, V, T):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6), sharex=True)
        # Gráfica de u(t)
        ax1.plot(T, U, color="blue")
        ax1.plot(T, U*0+uE, color="red")
        ax1.set_ylabel("u(t)")
        ax1.set_title("Evolución de u y v con Euler explícito")
        ax1.grid(True)

        # Gráfica de v(t)
        ax2.plot(T, V, color="blue")
        ax2.plot(T, V*0+vE, color="red")

        ax2.set_xlabel("Tiempo")
        ax2.set_ylabel("v(t)")
        ax2.grid(True)
        plt.tight_layout()
        plt.show()
    def grafDigramaEstabala():
        # Rango de valores
        a_vals = np.linspace(0.01, 2, 400)  # evitar división por cero
        b_vals = np.linspace(0.01, 2, 400)
        A, B = np.meshgrid(a_vals, b_vals)

        # Expresión común
        expr = 1 + (A+B)**2 - (2*B)/(A+B)

        # Discriminante
        disc = expr**2 - 4*(A+B)**2

        # Raíces lambda1, lambda2 (con seguridad para valores negativos en disc)
        sqrt_disc = np.sqrt(np.maximum(disc, 0))

        lambda1 = (-expr + sqrt_disc)/2
        lambda2 = (-expr - sqrt_disc)/2

        # Región donde AMBOS eigenvalores son negativos
        region = (lambda1 < 0) & (lambda2 < 0)

        # Graficar
        plt.figure(figsize=(6,5))
        plt.contourf(A, B, region, levels=[-0.5,0.5,1.5], colors=["white","lightblue"])
        plt.xlabel("a")
        plt.ylabel("b")
        #plt.title("Región donde λ1 < 0 y λ2 < 0")
        #plt.colorbar(label="1 si λ1,λ2 < 0; 0 si no")
        # Guardar antes de mostrar
        plt.savefig('GrafEstabilidadF.eps', format='eps', dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    def grafRelacionDispersion(p):
        eta = p[3]
        a = p[0]
        b = p[1]
        d= p[2]
        J, (ue, ve) = Jacobiano_en_equilibrio(p)
        fu=J[0,0]
        fv=J[0,1]
        gu=J[1,0]
        gv=J[1,1]
        c2 = np.linspace(-0.1, 20, 1000)

        # Cálculos
        Traza = fu + gv - (d + 1) * c2
        Det = d * c2**2 - c2 * (fu + gv * d) + fu * gv - fv * gu

        # Usar np.sqrt (en lugar de math.sqrt)
        Discriminante = Traza**2 - 4 * Det
        RelDis = (Traza + np.sqrt(Discriminante)) / 2  # signo +, rama de crecimiento

        # Gráfico
        plt.plot(c2, RelDis, label='d<d_c')


        d=0.05327625555898607
        Traza = fu + gv - (d + 1) * c2
        Det = d * c2**2 - c2 * (fu + gv * d) + fu * gv - fv * gu

        # Usar np.sqrt (en lugar de math.sqrt)
        Discriminante = Traza**2 - 4 * Det
        RelDis = (Traza + np.sqrt(Discriminante)) / 2  # signo +, rama de crecimiento

        # Gráfico
        plt.plot(c2, RelDis, label='d=d_c')

        d=0.05327625555898607+0.01
        Traza = fu + gv - (d + 1) * c2
        Det = d * c2**2 - c2 * (fu + gv * d) + fu * gv - fv * gu

        # Usar np.sqrt (en lugar de math.sqrt)
        Discriminante = Traza**2 - 4 * Det
        RelDis = (Traza + np.sqrt(Discriminante)) / 2  # signo +, rama de crecimiento

        # Gráfico
        plt.plot(c2, RelDis, label='d>d_c')
        plt.axhline(0, color='black', linestyle='--')
        plt.xlabel('Número de onda (c²)')
        plt.ylabel('Re(λ)')
        plt.legend()
        plt.grid(True)
        plt.savefig('RelacionDispersion.eps', format='eps', dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    def dcritico(self):
        J, (ue, ve) = self.Jacobiano_en_equilibrio(self.p)
        fu=J[0,0]
        fv=J[0,1]
        gu=J[1,0]
        gv=J[1,1]
        #print('jacobiano',J)
        h=2*fu*gv -4*(fu*gv-fv*gu)
        h=-2*fu*gv+4*fv*gu
        dc=(-h- math.sqrt(h**2-4*(fu**2)*(gv**2)))/(2*gv**2)
        return dc
    def etacritico(self,NOnda):
        J, (ue, ve) = self.Jacobiano_en_equilibrio(self.p)
        fu=J[0,0]
        fv=J[0,1]
        gu=J[1,0]
        gv=J[1,1]

        dc=self.dcritico()

        A = fu * gv - fv * gu
        B = -NOnda * (fu + gv * dc)
        C = dc * NOnda**2

        if abs(A) < 1e-12:
            etax1 = C / (-B)
            etax2 = C / (-B)

        else:
            disc = B**2 - 4*A*C

            if disc < 0:
                etax1 = np.nan
                etax2 = np.nan
            else:
                sqrt_disc = np.sqrt(disc)
                etax1 = (-B + sqrt_disc) / (2*A)
                etax2 = (-B - sqrt_disc) / (2*A)

        return etax1,etax2