# Turing Patterns in Curves

Numerical simulation of reaction--diffusion systems and Turing pattern formation on periodic curves and heterogeneous geometries.

---

## Description

This project implements a computational framework for studying:

- Turing instabilities
- Reaction--diffusion systems
- Heterogeneous diffusion
- Asymptotic homogenization
- Curvature-induced effects
- Pattern formation on parametrized curves

The code allows numerical simulations and visualization of patterns generated on periodic geometries embedded in $\mathbb{R}^3$.

---

## Project Structure

```text
.
├── main.py
├── geometria.py
├── SistemaRD_Schnakenberg.py
├── homogeneizacion.py
└── Resultados/
```

---

## Requirements

Install the required libraries:

```bash
pip install numpy matplotlib scipy
```

---

## Geometry Definition

The geometry is constructed using:

```python
geomToro = geometria.geometria_class(1, p_geometria)
```

with

```python
p_geometria = [R, r, p, q, Sfin]
```

### Parameters

| Parameter | Description |
|---|---|
| `R` | Major radius |
| `r` | Minor radius |
| `p` | Oscillation frequency |
| `q` | High-frequency parameter |
| `Sfin` | Final parameter interval |

Example:

```python
Sfin = 2*np.pi
r = 0.3
R = 1
eps = 1/7
p = 40
q = 5/eps
```

---

## Geometry Visualization

### Curve Plot

```python
geomToro.grafica_curva(0, Sfin)
```

### Curvature Plot

```python
geomToro.GraficaCurvaura(0, Sfin)
```

### Tangent Vector Plot

```python
geomToro.GraficaVTangente(0, Sfin)
```

---

## Schnakenberg Reaction--Diffusion Model

The implemented system corresponds to the Schnakenberg model:

\[
\begin{aligned}
u_t &= \nabla \cdot (D_u \nabla u) + a - u + u^2v, \\
v_t &= d \nabla \cdot (D_v \nabla v) + b - u^2v.
\end{aligned}
\]

Example parameters:

```python
a = 0.1
b = 1.5
d = 0.052
eta = 0.15
alpha = 0.01
```

The model includes heterogeneous diffusion induced by curvature.

---

## Diffusion Visualization

The diffusion coefficients can be visualized using:

```python
SisRD.GraficaDifusion(0, Sfin, num_points=10000)
```

---

## Homogenization

Effective diffusion coefficients are computed through asymptotic homogenization:

```python
promu, promv = H_RD.Coef_Efectivo(0, 2*np.pi)
```

This allows comparison between:

- Heterogeneous systems
- Homogenized effective systems

---

## Turing Instability Analysis

The code includes numerical exploration of Turing instability regions:

```python
H_RD.InestabilidadHomogeneizado(
    inter_d,
    inter_alpha,
    labelx="d",
    labely="alpha"
)
```

Parameters explored:

```python
d_min = 0.02
d_max = 0.06

alpha_min = 0.05
alpha_max = 1.5
```

---

## Numerical Simulations

The framework performs simulations of:

- Full heterogeneous models
- Homogenized models

### Heterogeneous Simulation

```python
H_RD.SimulRD(
    p,
    S=SisRD.geom.p_geom[-1],
    num_points=600,
    times_to_save=times_to_save,
    periodic=True
)
```

### Homogenized Simulation

```python
H_RD.SimulHomogenizacion(
    p,
    S=SisRD.geom.p_geom[-1],
    num_points=600,
    times_to_save=times_to_save,
    periodic=True
)
```

### Saved Times

```python
times_to_save = [
    0,10,20,30,40,50,
    100,150,200,
    500,600,700,800,900,1000
]
```

---

## Solution Comparison

The heterogeneous and homogenized solutions can be compared through:

```python
H_RD.GrafComparacion(
    FilenameU=NombreDelArchivoU,
    FilenameH=NombreDelArchivoH,
    FileSalida=Salida
)
```

---

## 3D Visualization of Solutions

The numerical solution can be projected over the geometry:

```python
geomToro.GrafColorF(
    Filename=NombreDelArchivo,
    FileSalida=NombreSalida,
    NombreFuncion="u(s,t)"
)
```

This produces 3D colored visualizations of the solution over the curve.

---

## Curvature Visualization in 3D

The curvature distribution can also be visualized:

```python
geomToro.GrafColorF(
    Filename="Curvatura.txt",
    FileSalida="Curvatura3D.png",
    NombreFuncion="kg(s)"
)
```

---

## Projection on the Plane

The solution can be projected onto a plane:

```python
geomToro.proyeccionSobreElPlano(
    Filename=NombreDelArchivo,
    FileSalida=NombreSalida,
    NombreFuncion=NombreFuncion0
)
```

---

## GIF Animations

Animations of pattern evolution can be generated:

```python
geomToro.generar_gif(
    1,
    "u(s,t)",
    filename="Resultados/solucionaHetero",
    salida_gif="Resultados/animacionHetero.gif",
    tiempos=times_to_save
)
```

---

## Output Files

All generated files are automatically stored in:

```text
Resultados/
```

including:

- Numerical solutions
- Curvature data
- Images
- GIF animations
- Comparison plots
- Copy of `main.py`

---

## Example Applications

This framework allows the study of:

- Turing pattern formation
- Curvature-induced diffusion effects
- Pattern localization
- Heterogeneous media
- Homogenization approximations
- Spatial oscillatory structures

---

## Future Work

Possible extensions include:

- Reaction--diffusion systems on surfaces
- Stability analysis on manifolds
- Spectral analysis of Laplace--Beltrami operators
- Hexagonal and stripe pattern classification
- Multi-species systems
- Nonlinear geometric effects

---

## Research project focused on:

- Reaction--diffusion systems
- Turing instability
- Geometric analysis
- Homogenization theory
- Pattern formation on periodic curves

## Authors

- Elkinn Adrian Calderon Barreto
- Dayana Carolina Colimba Tarapuez
