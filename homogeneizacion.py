import numpy as np
import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize._numdiff import approx_derivative
import math

class H_RD_Class:
    def __init__(self,sistema_RD):
        self.sistem=sistema_RD 
        self.hat_du,self.hat_dv=self.Coef_Efectivo(0,self.sistem.geom.p_geom[-1])

    def Coef_Efectivo(self,a,b,Param1=None,Param2=None,num_points=1000):
        intT = np.linspace(a, b, num_points)
        
        h=(b-a)/num_points

        bar_du_a,bar_dv_a=self.sistem.difusion(a,Param1,Param2)
        bar_du_b,bar_dv_b=self.sistem.difusion(b,Param1,Param2)
        
        dx,dy,dz,sig_a=self.sistem.geom.VectorTangente(a)
        dx,dy,dz,sig_b=self.sistem.geom.VectorTangente(b)

        bar_du_a=1.0/(bar_du_a/sig_a)
        bar_du_b=1.0/(bar_du_b/sig_b)

        bar_dv_a=1.0/(bar_dv_a/sig_a)
        bar_dv_b=1.0/(bar_dv_b/sig_b)


        hat_du=0.5*(bar_du_a+bar_du_b)
        hat_dv=0.5*(bar_dv_a+bar_dv_b)

        for i in range(1, len(intT)-1):
            t =a+ i*h
            
            bar_du_t,bar_dv_t=self.sistem.difusion(t,Param1,Param2)

            dx,dy,dz,sig_t=self.sistem.geom.VectorTangente(t)

            bar_du_t=1.0/(bar_du_t/sig_t)
            bar_dv_t=1.0/(bar_dv_t/sig_t)

            hat_du= hat_du+bar_du_t
            hat_dv= hat_dv+bar_dv_t
        
        

        hat_du=(h*hat_du)/((b-a))
        hat_dv=(h*hat_dv)/((b-a))

        hat_du=1.0/hat_du
        hat_dv=1.0/hat_dv

        sigma_du=self.sistem.geom.LongitudArco(a,b,num_points)[0]/(b-a)
        sigma_dv=self.sistem.geom.LongitudArco(a,b,num_points)[0]/(b-a) 

        hat_du=hat_du/sigma_du
        hat_dv=hat_dv/sigma_dv

        return hat_du,hat_dv
    
    def SimulHomogenizacion(self,p,S=2 * np.pi,num_points=1000,times_to_save=[0, 100],periodic=False,FileName="solucion",dte=1e-4):
        from pde import PDE, FieldCollection, ScalarField, CartesianGrid
        from pde import MemoryStorage
        import numpy as np
        import matplotlib.pyplot as plt
        from pde import ProgressTracker

        np.random.seed(123)

        print("Parámetros del sistema:", p)
        # -------------------------------------------------
        # malla espacial
        # -------------------------------------------------

        grid = CartesianGrid(
            [[0, S]],
            [num_points],
            periodic=periodic
        )

        x = grid.cell_coords[:, 0]

        # -------------------------------------------------
        # difusión heterogénea / homogenizada
        # -------------------------------------------------

        
        hatU, hatV = self.Coef_Efectivo(0, S, num_points=num_points)

        def Du_func(x):
            return hatU

        def Dv_func(x):
            return hatV


        D_u = ScalarField(
            grid,
            data=Du_func(x),
            label="D_u"
        )

        D_v = ScalarField(
            grid,
            data=Dv_func(x),
            label="D_v"
        )

        # -------------------------------------------------
        # PDE
        # -------------------------------------------------

        bc_type = None if periodic else "neumann"
        divergente=True;
        
        eq = PDE(
            {
                "u": "D_u*laplace(u) + eta*(a - u + u**2*v)",
                "v": "D_v*laplace(v) + eta*(b - u**2*v)",
            },
            consts={
                "a": p[0],
                "b": p[1],
                "eta": p[3],
                "D_u": D_u,
                "D_v": D_v,
            },
            bc=bc_type
        )

        # -------------------------------------------------
        # condición inicial
        # -------------------------------------------------

        Ue = self.sistem.Equilibrio(p)[0]
        Ve = self.sistem.Equilibrio(p)[1]
        print('equilibrios  H:', Ue, Ve)

        u = ScalarField(
            grid,
            data=Ue + 0.1 * np.random.normal(size=grid.shape),
            label="u"
        )

        v = ScalarField(
            grid,
            data=Ve + 0.08 * np.random.normal(size=grid.shape),
            label="v"
        )

        state = FieldCollection([u, v])

        # -------------------------------------------------
        # simulación y guardado
        # -------------------------------------------------

        storage = MemoryStorage()

        eq.solve(
            state,
            t_range=times_to_save[-1],
            dt=dte,
            tracker=[
                storage.tracker(interrupts=times_to_save),
                ProgressTracker()
            ]
        )
        # -------------------------------------------------
        # guardar resultados + gráficas
        # -------------------------------------------------

        for i, (t, fields) in enumerate(storage.items()):

            # ----------------------------------
            # guardar txt
            # ----------------------------------

            datos = np.column_stack((
                x,
                fields[0].data,
                fields[1].data
            ))

            nombre_archivo = f"{FileName}_{t:.2f}.txt"

            np.savetxt(
                nombre_archivo,
                datos,
                header="x u v",
                fmt="%.8f"
            )

            # ----------------------------------
            # gráficas
            # ----------------------------------

            fig, (ax1, ax2) = plt.subplots(
                2, 1,
                figsize=(8, 6),
                sharex=True
            )

            # u
            ax1.plot(x, fields[0].data)
            ax1.axhline(Ue, linestyle="--")
            ax1.set_ylabel("u(x,t)")
            ax1.set_title(f"Solución en t = {t:.2f}")
            ax1.grid(True)

            # v
            ax2.plot(x, fields[1].data)
            ax2.axhline(Ve, linestyle="--")
            ax2.set_xlabel("x")
            ax2.set_ylabel("v(x,t)")
            ax2.grid(True)

            plt.tight_layout()
            plt.show()

        return storage

    def SimulRD(self,p,S=2 * np.pi,num_points=1000,times_to_save=[0, 100],periodic=False,FileName="solucion",dte=1e-4):
        from pde import PDE, FieldCollection, ScalarField, CartesianGrid
        from pde import MemoryStorage
        import numpy as np
        import matplotlib.pyplot as plt
        from pde import ProgressTracker

        np.random.seed(123)

        print("Parámetros del sistema:", p)
        # -------------------------------------------------
        # malla espacial
        # -------------------------------------------------

        grid = CartesianGrid(
            [[0, S]],
            [num_points],
            periodic=periodic
        )

        x = grid.cell_coords[:, 0]

        # -------------------------------------------------
        # difusión heterogénea / homogenizada
        # -------------------------------------------------


        def Du_func(x):
            return self.sistem.difusion(x)[0]

        def Dv_func(x):
            return self.sistem.difusion(x)[1]

        def SigmaU_func(x):
            return self.sistem.geom.VectorTangente(x)[3]

        def SigmaV_func(x):
            return self.sistem.geom.VectorTangente(x)[3]

        D_u = ScalarField(
            grid,
            data=Du_func(x),
            label="D_u"
        )

        D_v = ScalarField(
            grid,
            data=Dv_func(x),
            label="D_v"
        )

        sigmaU = ScalarField(
            grid,
            data=SigmaU_func(x),
            label="sigmaU"
        )

        sigmaV = ScalarField(
            grid,
            data=SigmaV_func(x),
            label="sigmaV"
        )

        # -------------------------------------------------
        # PDE
        # -------------------------------------------------

        bc_type = None if periodic else "neumann"
        divergente=True;
        if divergente:
            eq = PDE(
                {
                    "u": "(1.0/sigmaU)*divergence((D_u/sigmaU) * gradient(u)) + eta*(a - u + u**2*v)",
                    "v": "(1.0/sigmaV)*divergence((D_v/sigmaV) * gradient(v)) + eta*(b - u**2*v)",
                },
                consts={
                    "a": p[0],
                    "b": p[1],
                    "eta": p[3],
                    "sigmaU": sigmaU,
                    "sigmaV": sigmaV,
                    "D_u": D_u,
                    "D_v": D_v,
                },
                bc=bc_type
            )
        else:

            eq = PDE(
                {
                    "u": "(1.0/sigmaU**2) * (D_u * laplace(u) +(1.0/sigmaU) *  dot(gradient(D_u/sigmaU), gradient(u))) + eta * (a - u + u**2 * v)",

                    "v": "(1.0/sigmaV**2) * (D_v * laplace(v) +(1.0/sigmaV) *  dot(gradient(D_v/sigmaU), gradient(v))) + eta * (b - u**2 * v)",
                },
                consts={
                    "a": p[0],
                    "b": p[1],
                    "eta": p[3],
                    "sigmaU": sigmaU,
                    "sigmaV": sigmaV,
                    "D_u": D_u,
                    "D_v": D_v,
                },
                bc=bc_type
            )

        # -------------------------------------------------
        # condición inicial
        # -------------------------------------------------

        Ue = self.sistem.Equilibrio(p)[0]
        Ve = self.sistem.Equilibrio(p)[1]
        print('equilibrios  H:', Ue, Ve)

        u = ScalarField(
            grid,
            data=Ue + 0.1 * np.random.normal(size=grid.shape),
            label="u"
        )

        v = ScalarField(
            grid,
            data=Ve + 0.08 * np.random.normal(size=grid.shape),
            label="v"
        )

        state = FieldCollection([u, v])

        # -------------------------------------------------
        # simulación y guardado
        # -------------------------------------------------

        storage = MemoryStorage()

        eq.solve(
            state,
            t_range=times_to_save[-1],
            dt=dte,
            tracker=[
                storage.tracker(interrupts=times_to_save),
                ProgressTracker()
            ]
        )
        # -------------------------------------------------
        # guardar resultados + gráficas
        # -------------------------------------------------

        for i, (t, fields) in enumerate(storage.items()):

            # ----------------------------------
            # guardar txt
            # ----------------------------------

            datos = np.column_stack((
                x,
                fields[0].data,
                fields[1].data
            ))

            nombre_archivo = f"{FileName}_{t:.2f}.txt"

            np.savetxt(
                nombre_archivo,
                datos,
                header="x u v",
                fmt="%.8f"
            )

            # ----------------------------------
            # gráficas
            # ----------------------------------

            fig, (ax1, ax2) = plt.subplots(
                2, 1,
                figsize=(8, 6),
                sharex=True
            )

            # u
            ax1.plot(x, fields[0].data)
            ax1.axhline(Ue, linestyle="--")
            ax1.set_ylabel("u(x,t)")
            ax1.set_title(f"Solución en t = {t:.2f}")
            ax1.grid(True)

            # v
            ax2.plot(x, fields[1].data)
            ax2.axhline(Ve, linestyle="--")
            ax2.set_xlabel("x")
            ax2.set_ylabel("v(x,t)")
            ax2.grid(True)

            plt.tight_layout()
            plt.show()

        return storage

    def GrafComparacion(self,FilenameU=None,FilenameH=None,FileSalida='comparacion.png'):
        
        if FilenameU is not None and FilenameH is not None: 
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            dataU = np.loadtxt(FilenameU)
            dataH = np.loadtxt(FilenameH)

            # columnas
            Ut = dataU[:, 0]
            Uu = dataU[:, 1]
            Uv = dataU[:, 2]

            # columnas
            Ht  = dataH[:, 0]
            Hu = dataH[:, 1]
            Hv = dataH[:, 2]


            # ----------------------------------
            # gráficas
            # ----------------------------------

            fig, (ax1, ax2) = plt.subplots(
                2, 1,
                figsize=(8, 6),
                sharex=True
            )
            Ue = self.sistem.Equilibrio(self.sistem.p)[0]
            Ve = self.sistem.Equilibrio(self.sistem.p)[1]
            # u
            ax1.plot(Ut, Uu)
            ax1.plot(Ht, Hu)
            ax1.set_ylabel("u(s,t)")
            ax1.grid(True)

            # v
            ax2.plot(Ut, Uv)
            ax2.plot(Ht, Hv)
            ax2.set_xlabel("s")
            ax2.set_ylabel("v(s,t)")
            ax2.grid(True)

            plt.tight_layout()
            plt.savefig(FileSalida, dpi=300, bbox_inches="tight")
            plt.show()
            
        else:
            print("Ingrese los dos archivos a comparar")
        



    def InestabilidadHomogeneizado(self,IntPara1,IntPara2,labelx='p1',labely='p2',FileSalida='InestabilidadTuring.png'):
        
        dc=self.sistem.dcritico()
        etac=self.sistem.etacritico(2*np.pi/self.sistem.geom.p_geom[-1])
        print('Difusion Critica:',dc)
        print('Etacritico para numero de onda 1:',etac)

        CndTuring = np.zeros((len(IntPara2), len(IntPara1)))

        for i, ip1 in enumerate(IntPara1):
            for j, ip2 in enumerate(IntPara2):

                hat_du, hat_dv = self.Coef_Efectivo(
                    0,
                    self.sistem.geom.p_geom[-1],
                    ip1,
                    ip2
                )

                if dc > hat_du / hat_dv:
                    CndTuring[j, i] = 1
                else:
                    CndTuring[j, i] = 0
        np.save("CndTuring.npy", CndTuring)

        fig, ax = plt.subplots(figsize=(8, 6))

        im = ax.contourf(
            IntPara1,
            IntPara2,
            CndTuring,
            levels=100,        # más niveles = más suavidad
            cmap="viridis"
        )

        ax.set_xlabel(labelx)
        ax.set_ylabel(labely)
        ax.set_title("Región de inestabilidad de Turing")

        plt.colorbar(im, ax=ax, label="CndTuring")
        plt.savefig(FileSalida, dpi=300, bbox_inches="tight")
        plt.show()



