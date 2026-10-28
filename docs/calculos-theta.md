# Obtención Paso a Paso de los Elementos de la Matriz para el Cálculo de los Parámetros $\Theta_{mn}$ (Método SGF 1D)

Este documento detalla el procedimiento matemático continuo y matricial para calcular los parámetros auxiliares $\Theta_{mn}$ del método **Spectral Green's Function (SGF)** en geometría 1D con ordenadas discretas ($S_N$), partiendo de los autovalores espectrales ($\nu_k$) y las autofunciones angulares ($a_n(\nu_k)$) ya conocidos.

---

## 1. Fundamento y Ecuación Auxiliar del Método SGF

El método SGF busca preservar exactamente la solución analítica espectral dentro de cada nodo espacial $i \in [x_{i-1/2}, x_{i+1/2}]$ de ancho $h_i = x_{i+1/2} - x_{i-1/2}$.

La **solución analítica general** (expresión espectral) para la dirección angular $m$ viene dada por:

$$\Psi_m(x) = \frac{Q}{2\sigma_t(1 - c_0)} + \sum_{k=1}^N \alpha_k a_m(\nu_k) e^{\frac{\sigma_t x}{\nu_k}}$$

donde:
* $Q$: fuente externa isotrópica fija.
* $\sigma_t$: sección eficaz total del nodo.
* $c_0 = \frac{\sigma_{s0}}{\sigma_t}$: ratio de dispersión (*scattering ratio*).
* $\nu_k$: autovalores del espectro espacial ($k = 1, \dots, N$).
* $a_m(\nu_k)$: componentes del autovector angular asociadas a la dirección $m$ y al modo $k$.
* $\alpha_k$: coeficientes determinados por las condiciones de contorno globales.

La **ecuación auxiliar del método SGF** relaciona el flujo angular promedio en el nodo $\Psi_{m,i}$ con los flujos angulares incidentes en las aristas/fronteras del nodo:

$$\Psi_{m,i} = \sum_{\mu_n > 0} \Theta_{mn} \Psi_{n,i-1/2} + \sum_{\mu_n < 0} \Theta_{mn} \Psi_{n,i+1/2} + G_m(Q)$$

---

## 2. Definición de Flujos y Sustitución de la Solución Espectral

### 2.1. Flujo Angular Promedio en el Nodo ($\Psi_{m,i}$)
Por definición de valor promedio espacial en el intervalo $[x_{i-1/2}, x_{i+1/2}]$:

$$\Psi_{m,i} = \frac{1}{h_i} \int_{x_{i-1/2}}^{x_{i+1/2}} \Psi_m(x) \, dx$$

Sustituyendo la solución analítica general en la integral:

$$\Psi_{m,i} = \frac{Q}{2\sigma_t(1 - c_0)} + \frac{1}{h_i} \sum_{k=1}^N \alpha_k a_m(\nu_k) \int_{x_{i-1/2}}^{x_{i+1/2}} e^{\frac{\sigma_t x}{\nu_k}} \, dx$$

Evaluando la integral exponencial:

$$\Psi_{m,i} = \frac{Q}{2\sigma_t(1 - c_0)} + \frac{1}{h_i} \sum_{k=1}^N \alpha_k a_m(\nu_k) \frac{\nu_k}{\sigma_t} \left[ e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} - e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \right] \quad \text{--- (Eq. 1)}$$

### 2.2. Flujos Angulares en las Aristas ($\Psi_m(x_{i \pm 1/2})$)
Evaluando la solución espectral directamente en las posiciones de borde $x = x_{i+1/2}$ y $x = x_{i-1/2}$:

$$\Psi_n(x_{i-1/2}) = \frac{Q}{2\sigma_t(1 - c_0)} + \sum_{k=1}^N \alpha_k a_n(\nu_k) e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \quad \text{--- (Eq. 2)}$$

$$\Psi_n(x_{i+1/2}) = \frac{Q}{2\sigma_t(1 - c_0)} + \sum_{k=1}^N \alpha_k a_n(\nu_k) e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} \quad \text{--- (Eq. 3)}$$

---

## 3. Deducción de la Ecuación Modal de Balance

Sustituyendo **Eq. 1**, **Eq. 2** y **Eq. 3** dentro de la ecuación auxiliar SGF:

$$\frac{Q}{2\sigma_t(1-c_0)} + \frac{1}{h_i} \sum_{k=1}^N \alpha_k \frac{\nu_k}{\sigma_t} a_m(\nu_k) \left( e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} - e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \right) =$$
$$\sum_{\mu_n > 0} \Theta_{mn} \left[ \frac{Q}{2\sigma_t(1-c_0)} + \sum_{k=1}^N \alpha_k a_n(\nu_k) e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \right] + \sum_{\mu_n < 0} \Theta_{mn} \left[ \frac{Q}{2\sigma_t(1-c_0)} + \sum_{k=1}^N \alpha_k a_n(\nu_k) e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} \right] + G_m(Q)$$

### 3.1. Separación de la Solución Particular ($G_m(Q)$)
Agrupando los términos independientes de $\alpha_k$ (correspondientes a la fuente $Q$):

$$G_m(Q) = \frac{Q}{2\sigma_t(1 - c_0)} \left( 1 - \sum_{n=1}^N \Theta_{mn} \right)$$

### 3.2. Separación de la Solución Homogénea (Modos $k$)
Para la parte homogénea, igualamos los coeficientes asociados a cada constante arbitraria $\alpha_k$. Puesto que la ecuación debe cumplirse para cualquier $\alpha_k$, se debe satisfacer de forma independiente para cada autovalor $k = 1, \dots, N$:

$$a_m(\nu_k) \frac{\nu_k}{h_i \sigma_t} \left( e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} - e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \right) = e^{\frac{\sigma_t x_{i-1/2}}{\nu_k}} \sum_{\mu_n > 0} \Theta_{mn} a_n(\nu_k) + e^{\frac{\sigma_t x_{i+1/2}}{\nu_k}} \sum_{\mu_n < 0} \Theta_{mn} a_n(\nu_k)$$

---

## 4. Centrado en el Nodo e Introducción de Funciones Hiperbólicas

Para simplificar los exponentes exponenciales, trasladamos el origen del sistema de coordenadas al centro del nodo espacial $x_i$:

$$x_i = \frac{x_{i+1/2} + x_{i-1/2}}{2}$$

por lo cual:
$$x_{i+1/2} - x_i = \frac{h_i}{2} \quad \text{y} \quad x_{i-1/2} - x_i = -\frac{h_i}{2}$$

Multiplicamos toda la ecuación anterior por el factor $e^{-\frac{\sigma_t x_i}{\nu_k}}$:

$$a_m(\nu_k) \frac{\nu_k}{h_i \sigma_t} \left( e^{\frac{\sigma_t h_i}{2\nu_k}} - e^{-\frac{\sigma_t h_i}{2\nu_k}} \right) = e^{-\frac{\sigma_t h_i}{2\nu_k}} \sum_{\mu_n > 0} \Theta_{mn} a_n(\nu_k) + e^{\frac{\sigma_t h_i}{2\nu_k}} \sum_{\mu_n < 0} \Theta_{mn} a_n(\nu_k)$$

Utilizando la identidad de la función seno hiperbólico $\sinh(a) = \frac{e^a - e^{-a}}{2}$, el término de la izquierda se reduce a:

$$e^{\frac{\sigma_t h_i}{2\nu_k}} - e^{-\frac{\sigma_t h_i}{2\nu_k}} = 2 \sinh\left( \frac{\sigma_t h_i}{2\nu_k} \right)$$

Llegando a la **Ecuación Fundamental de Acoplamiento Modal (Eq. 8)**:

$$\mathbf{a_m(\nu_k) \frac{2\nu_k}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_k} \right) = e^{-\frac{\sigma_t h_i}{2\nu_k}} \sum_{\mu_n > 0} \Theta_{mn} a_n(\nu_k) + e^{\frac{\sigma_t h_i}{2\nu_k}} \sum_{\mu_n < 0} \Theta_{mn} a_n(\nu_k)} \quad \text{--- (Eq. 4)}$$

---

## 5. Construcción del Sistema Matricial de Ecuaciones

Fijando una dirección angular $m \in \{1, \dots, N\}$, la **Eq. 4** debe evaluarse para cada autovalor $k = 1, \dots, N$. Esto genera un sistema de $N$ ecuaciones lineales con $N$ incógnitas, donde las incógnitas son la fila de parámetros $\Theta_{m1}, \Theta_{m2}, \dots, \Theta_{mN}$.

### 5.1. Regla de Signos para las Columnas de la Matriz $\mathbf{M}$
Para cada dirección discreta $n$ ($n = 1, \dots, N$):
* Si la dirección es **hacia la derecha** ($\mu_n > 0$), el coeficiente asociado en la fila $k$ es:
  $$M_{k,n} = e^{-\frac{\sigma_t h_i}{2\nu_k}} a_n(\nu_k)$$
* Si la dirección es **hacia la izquierda** ($\mu_n < 0$), el coeficiente asociado en la fila $k$ es:
  $$M_{k,n} = e^{\frac{\sigma_t h_i}{2\nu_k}} a_n(\nu_k)$$

> **Nota Clave:** La matriz de coeficientes $\mathbf{M}$ de dimensión $N \times N$ es **idéntica para todas las direcciones $m$**. Lo único que cambia de una dirección $m$ a otra es el vector independiente del lado derecho $\mathbf{b}_m$.

### 5.2. Vector del Lado Derecho $\mathbf{b}_m$
Para una dirección fija $m$, la componente $k$-ésima del vector $\mathbf{b}_m$ es:

$$b_{m,k} = a_m(\nu_k) \frac{2\nu_k}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_k} \right)$$

### 5.3. Formulación Matricial General
El sistema lineal se expresa de forma compacta como:

$$\mathbf{M} \cdot \boldsymbol{\Theta}_m^T = \mathbf{b}_m$$

o explícitamente:

$$\begin{bmatrix}
M_{1,1} & M_{1,2} & \dots & M_{1,N} \\
M_{2,1} & M_{2,2} & \dots & M_{2,N} \\
\vdots & \vdots & \ddots & \vdots \\
M_{N,1} & M_{N,2} & \dots & M_{N,N}
\end{bmatrix}
\begin{bmatrix}
\Theta_{m1} \\
\Theta_{m2} \\
\vdots \\
\Theta_{mN}
\end{bmatrix}
=
\begin{bmatrix}
a_m(\nu_1) \frac{2\nu_1}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_1} \right) \\
a_m(\nu_2) \frac{2\nu_2}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_2} \right) \\
\vdots \\
a_m(\nu_N) \frac{2\nu_N}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_N} \right)
\end{bmatrix}$$

---

## 6. Caso Práctico Detallado: Cuadratura $S_2$ ($N=2$)

Para $N=2$, tenemos dos direcciones angulares:
* $\mu_1 < 0$ (izquierda)
* $\mu_2 > 0$ (derecha)

Y dos autovalores espectrales: $\nu_1$ y $\nu_2$.

### 6.1. Matriz $\mathbf{M}$ para $S_2$
Como $\mu_1 < 0$ (columna 1) y $\mu_2 > 0$ (columna 2):

$$\mathbf{M} = \begin{bmatrix}
e^{\frac{\sigma_t h_i}{2\nu_1}} a_1(\nu_1) & e^{-\frac{\sigma_t h_i}{2\nu_1}} a_2(\nu_1) \\[6pt]
e^{\frac{\sigma_t h_i}{2\nu_2}} a_1(\nu_2) & e^{-\frac{\sigma_t h_i}{2\nu_2}} a_2(\nu_2)
\end{bmatrix}$$

### 6.2. Sistema para calcular la primera fila de $\Theta$ ($m=1$: dirección $\mu_1 < 0$)

$$\begin{bmatrix}
e^{\frac{\sigma_t h_i}{2\nu_1}} a_1(\nu_1) & e^{-\frac{\sigma_t h_i}{2\nu_1}} a_2(\nu_1) \\[6pt]
e^{\frac{\sigma_t h_i}{2\nu_2}} a_1(\nu_2) & e^{-\frac{\sigma_t h_i}{2\nu_2}} a_2(\nu_2)
\end{bmatrix}
\begin{bmatrix}
\Theta_{11} \\[6pt]
\Theta_{12}
\end{bmatrix}
=
\begin{bmatrix}
a_1(\nu_1) \frac{2\nu_1}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_1} \right) \\[6pt]
a_1(\nu_2) \frac{2\nu_2}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_2} \right)
\end{bmatrix}$$

### 6.3. Sistema para calcular la segunda fila de $\Theta$ ($m=2$: dirección $\mu_2 > 0$)

$$\begin{bmatrix}
e^{\frac{\sigma_t h_i}{2\nu_1}} a_1(\nu_1) & e^{-\frac{\sigma_t h_i}{2\nu_1}} a_2(\nu_1) \\[6pt]
e^{\frac{\sigma_t h_i}{2\nu_2}} a_1(\nu_2) & e^{-\frac{\sigma_t h_i}{2\nu_2}} a_2(\nu_2)
\end{bmatrix}
\begin{bmatrix}
\Theta_{21} \\[6pt]
\Theta_{22}
\end{bmatrix}
=
\begin{bmatrix}
a_2(\nu_1) \frac{2\nu_1}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_1} \right) \\[6pt]
a_2(\nu_2) \frac{2\nu_2}{h_i \sigma_t} \sinh\left( \frac{\sigma_t h_i}{2\nu_2} \right)
\end{bmatrix}$$

---

## 7. Algoritmo de Implementación Numérica (Pseudocódigo)

```python
import numpy as np

def calcular_matriz_theta(N, sigma_t, h_i, nu, a):
    """
    N: orden de la quadratura (ej. 2, 4, 8)
    sigma_t: sección eficaz total
    h_i: ancho del nodo espacial
    nu: array de autovalores de tamaño N [nu_1, nu_2, ..., nu_N]
    a: matriz de autofunciones de tamaño (N_direcciones, N_modos)
    """
    M = np.zeros((N, N))
    
    # 1. Construir la matriz de coeficientes M
    for k in range(N): # iterar sobre autovalores
        for n in range(N): # iterar sobre direcciones
            arg = (sigma_t * h_i) / (2.0 * nu[k])
            if mu[n] < 0:
                M[k, n] = np.exp(arg) * a[n, k]
            else: # mu[n] > 0
                M[k, n] = np.exp(-arg) * a[n, k]
                
    Theta = np.zeros((N, N))
    
    # 2. Resolver el sistema M * Theta_m = b_m para cada dirección m
    for m in range(N):
        b_m = np.zeros(N)
        for k in range(N):
            arg = (sigma_t * h_i) / (2.0 * nu[k])
            b_m[k] = a[m, k] * (2.0 * nu[k] / (h_i * sigma_t)) * np.sinh(arg)
            
        # Resolver sistema lineal M * Theta[m, :] = b_m
        Theta[m, :] = np.linalg.solve(M, b_m)
        
    return Theta
```

---

## 8. Resumen Final

1. La matriz de coeficientes $\mathbf{M}$ se ensambla evaluando la exponencial espacial centrada $e^{\pm \frac{\sigma_t h_i}{2\nu_k}}$ multiplicada por las autofunciones angulares $a_n(\nu_k)$.
2. El signo del exponente depende únicamente de la dirección espacial del flujo entrante en las caras del nodo ($\mu_n > 0$ o $\mu_n < 0$).
3. La resolución del sistema $\mathbf{M} \boldsymbol{\Theta}_m^T = \mathbf{b}_m$ proporciona directamente las constantes $\Theta_{mn}$ para el cálculo del flujo promedio $\Psi_{m,i}$ libre de errores de truncamiento de malla.
