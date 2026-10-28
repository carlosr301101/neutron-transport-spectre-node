# Resolución Analítica del Problema de Transporte Neutrónico 1D ($S_4$) mediante Análisis Espectral

## 1. Definición y Datos del Problema

Se considera el mismo medio homogéneo unidimensional (placa plana) del ejercicio $S_2$: sin fuentes internas fijas ($Q = 0$), sin fisión, con flujo neutrónico incidente por la frontera izquierda. La única diferencia es el **aumento del orden de la cuadratura angular**, de $S_2$ ($N = 2$) a $S_4$ ($N = 4$), lo que duplica el número de direcciones y, por tanto, el número de modos espectrales.

### Parámetros Físicos y Geométricos
* **Espesor de la placa ($L$):** $100\text{ cm}$ ($x \in [0, 100]$)
* **Sección eficaz total ($\sigma_t$):** $1.0\text{ cm}^{-1}$
* **Sección eficaz de dispersión isotrópica ($\sigma_{s0}$):** $0.97\text{ cm}^{-1}$
* **Razón de dispersión / Albedo de colisión ($c$):** $c = \frac{\sigma_{s0}}{\sigma_t} = \frac{0.97}{1.0} = 0.97$

### Condiciones de Contorno
* **Frontera Izquierda ($x = 0\text{ cm}$):** Flujo entrante conocido $\psi_m(0) = 1.0$ para todas las direcciones con $\mu_m > 0$ (es decir, $m = 3, 4$).
* **Frontera Derecha ($x = 100\text{ cm}$):** Condición de vacío $\psi_m(100) = 0.0$ para todas las direcciones con $\mu_m < 0$ (es decir, $m = 1, 2$).

En $S_4$ hay **dos** direcciones entrantes y **dos** salientes en cada frontera, por lo que se imponen cuatro condiciones de contorno para determinar los cuatro modos de la solución.

---

## 2. Discretización Angular en Ordenadas Discretas ($S_4$)

Para la ordenación discreta $S_4$ ($N = 4$), los puntos de evaluación $\mu_m$ y los pesos $w_m$ provienen de la cuadratura de Gauss-Legendre (raíces del polinomio $P_4(\mu) = 0$), cuyos valores se toman del archivo `cuadraturas.py`:

$$\{\mu_k\}_{\mu>0} = \{\,0.861136311594053,\; 0.339981043584856\,\}$$
$$\{w_k\}_{\mu>0} = \{\,0.347854845137454,\; 0.652145154862546\,\}$$

En forma cerrada exacta, estas son:

$$\mu = \pm \sqrt{\frac{3 \mp 2\sqrt{6/5}}{7}}, \qquad w = \frac{18 \pm \sqrt{30}}{36}$$

Como en $S_2$, el conjunto de direcciones es simétrico ($\mu_{m'} = -\mu_m$ para el índice espejo) y los pesos son iguales para cada par $\pm\mu_k$. Reordenando los índices de izquierda a derecha, se tiene:

| $m$ | $\mu_m$ | $w_m$ | Sentido |
|:---:|:---:|:---:|:---|
| 1 | $-0.8611363$ | $0.3478548$ | hacia la izquierda |
| 2 | $-0.3399810$ | $0.6521452$ | hacia la izquierda |
| 3 | $+0.3399810$ | $0.6521452$ | hacia la derecha |
| 4 | $+0.8611363$ | $0.3478548$ | hacia la derecha |

* **Pesos de integración:** $w_1 + w_2 + w_3 + w_4 = 2.0$ (igual que en $S_2$).

El sistema acoplado de cuatro ecuaciones diferenciales ordinarias $S_4$ es:

$$\mu_m \frac{d\psi_m(x)}{dx} + \sigma_t \psi_m(x) = \frac{\sigma_{s0}}{2} \sum_{n=1}^{4} w_n \psi_n(x), \quad m = 1, 2, 3, 4$$

es decir, explícitamente:

$$\mu_1 \frac{d\psi_1}{dx} + \sigma_t \psi_1 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_2 \frac{d\psi_2}{dx} + \sigma_t \psi_2 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_3 \frac{d\psi_3}{dx} + \sigma_t \psi_3 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_4 \frac{d\psi_4}{dx} + \sigma_t \psi_4 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$

donde el flujo escalar asociado a esta cuadratura es

$$\phi(x) = \sum_{m=1}^{4} w_m \psi_m(x), \qquad \sum_{m=1}^{4} w_m = 2$$

Este es el factor que aparece en el término de fuente $\frac{\sigma_{s0}}{2}\phi(x)$; como en $S_2$, no hay un factor $1/2$ adicional.

---

## 3. Separación de Variables (Ansatz Espectral Modal)

Planteamos la solución mediante modos exponenciales homogéneos de la forma:

$$\psi_m(x) = v_m(\nu) e^{-x/\nu}, \quad m = 1, 2, 3, 4$$

donde:
* $\nu$ es el autovalor espectral (longitud de relajación o atenuación espacial).
* $v_m(\nu)$ es la componente $m$-ésima del autovector angular correspondiente.

Sustituyendo el ansatz en el sistema de ecuaciones y usando que $\frac{d\psi_m}{dx} = -\frac{1}{\nu} v_m e^{-x/\nu}$:

$$\left( \sigma_t - \frac{\mu_m}{\nu} \right) v_m(\nu) = \frac{\sigma_{s0}}{2} \underbrace{\sum_{n=1}^{4} w_n v_n(\nu)}_{\text{flujo escalar modal}}$$

Normalizando el flujo escalar modal a la unidad ($\sum_n w_n v_n(\nu) = 1$) y recordando que $\sigma_{s0}/\sigma_t = c$, se despeja el autovector angular:

$$v_m(\nu) = \frac{c/2}{1 - \frac{\mu_m}{\sigma_t \nu}} = \frac{0.485}{1 - \frac{\mu_m}{\nu}} \qquad (\sigma_t = 1.0\text{ cm}^{-1})$$

A diferencia de $S_2$ (un solo par de autovalores $\pm\nu$), en $S_4$ la condición de consistencia conduce a una ecuación de **cuarto grado**, por lo que existirán **dos pares** de autovalores $\pm\nu_1$ y $\pm\nu_2$ y, con ellos, **cuatro modos** espectrales.

---

## 4. Deducción Analítica de los Autovalores ($\nu$)

### 4.1 Condición de consistencia

Imponiendo $\sum_n w_n v_n(\nu) = 1$ en la fórmula del autovector:

$$\frac{c}{2} \sum_{m=1}^{4} \frac{w_m}{1 - \frac{\mu_m}{\sigma_t \nu}} = 1 \quad \Longrightarrow \quad \sum_{m=1}^{4} \frac{w_m}{1 - \frac{\mu_m}{\sigma_t \nu}} = \frac{2}{c}$$

### 4.2 Reducción por simetría

Como $\mu_{1,2} = -\mu_{3,4}$ con pesos idénticos, se agrupan los pares $\pm\mu_k$ ($k=1,2$):

$$\sum_{k=1}^{2} w_k \left[ \frac{1}{1 - \frac{\mu_k}{\sigma_t \nu}} + \frac{1}{1 + \frac{\mu_k}{\sigma_t \nu}} \right] = \frac{2}{c}$$

La suma entre corchetes es $\dfrac{2}{1 - \left(\frac{\mu_k}{\sigma_t \nu}\right)^2}$, de modo que

$$\sum_{k=1}^{2} \frac{w_k}{1 - \frac{\mu_k^2}{(\sigma_t \nu)^2}} = \frac{1}{c}$$

y con $\sigma_t = 1.0\text{ cm}^{-1}$ (numéricamente $\sigma_t\nu = \nu$):

$$\frac{w_1 \nu^2}{\nu^2 - \mu_1^2} + \frac{w_2 \nu^2}{\nu^2 - \mu_2^2} = \frac{1}{c}$$

### 4.3 Ecuación característica

Multiplicando por $(\nu^2 - \mu_1^2)(\nu^2 - \mu_2^2)$ y usando $w_1 + w_2 = 1$:

$$w_1 \nu^2 (\nu^2 - \mu_2^2) + w_2 \nu^2 (\nu^2 - \mu_1^2) = \frac{1}{c} (\nu^2 - \mu_1^2)(\nu^2 - \mu_2^2)$$

$$\nu^4 - Q \nu^2 = \frac{1}{c} \Big( \nu^4 - P \nu^2 + R \Big)$$

donde se definieron las constantes (que dependen solo de la cuadratura):

$$P = \mu_1^2 + \mu_2^2, \qquad Q = w_1\mu_2^2 + w_2\mu_1^2, \qquad R = \mu_1^2\mu_2^2$$

Reagrupando se obtiene una **ecuación bicuadrática** en $\nu$:

$$(1 - c)\nu^4 + (cQ - P)\nu^2 + R = 0$$

### 4.4 Valores numéricos exactos

Con $\mu_1^2 = 0.7415557$, $\mu_2^2 = 0.1155871$, $w_1 = 0.3478548$, $w_2 = 0.6521452$:

$$P = \mu_1^2 + \mu_2^2 = 0.8571429 = \frac{6}{7}$$

$$Q = w_1\mu_2^2 + w_2\mu_1^2 = 0.5238095 = \frac{11}{21}$$

$$R = \mu_1^2\mu_2^2 = 0.0857143 = \frac{3}{35}$$

$$cQ - P = 0.97 \cdot \frac{11}{21} - \frac{6}{7} = \frac{10.67}{21} - \frac{18}{21} = -\frac{733}{2100} \approx -0.3490476$$

Sustituyendo en la ecuación característica:

$$\frac{3}{100}\nu^4 - \frac{733}{2100}\nu^2 + \frac{3}{35} = 0$$

Multiplicando todo por $2100$ se obtiene la **ecuación característica exacta de $S_4$**:

$$\boxed{\;63\nu^4 - 733\nu^2 + 180 = 0\;}$$

(Compárese con $S_2$, donde la condición se reducía a la simple $\nu^2 = 100/9$.)

### 4.5 Raíces

Con $y = \nu^2$, la ecuación se reduce a la cuadrática $63y^2 - 733y + 180 = 0$:

$$y = \frac{733 \pm \sqrt{733^2 - 4 \cdot 63 \cdot 180}}{2 \cdot 63} = \frac{733 \pm \sqrt{491929}}{126}$$

$$\sqrt{491929} = 701.3765037$$

$$y_1 = \nu_1^2 = \frac{733 + 701.3765037}{126} = 11.3839405 \quad \Longrightarrow \quad \nu_1 = 3.3740096\text{ cm}$$
$$y_2 = \nu_2^2 = \frac{733 - 701.3765037}{126} = 0.2509801 \quad \Longrightarrow \quad \nu_2 = 0.5009792\text{ cm}$$

Obtenemos así **cuatro** autovalores simétricos (dos pares):

| Autovalor | $\nu_j$ | $1/\nu_j$ ($\text{cm}^{-1}$) | Tipo de modo |
|:---:|:---:|:---:|:---|
| $\nu_1$ | $+3.3740096$ | $0.2963833$ | fundamental (decae hacia la derecha) |
| $-\nu_1$ | $-3.3740096$ | $-0.2963833$ | fundamental (crece hacia la derecha) |
| $\nu_2$ | $+0.5009792$ | $1.9960910$ | transitorio (capa límite izquierda) |
| $-\nu_2$ | $-0.5009792$ | $-1.9960910$ | transitorio (capa límite derecha) |

**Observaciones:**
1. El autovalor fundamental $\nu_1 = 3.3740096\text{ cm}$ es mayor que el de $S_2$ ($10/3 = 3.3333\text{ cm}$): al refinar la cuadratura, la longitud de relajación aumenta, es decir, el flujo decae más lentamente.
2. El segundo autovalor $\nu_2 = 0.5009792\text{ cm}$ es mucho menor que el espesor: define un modo transitorio que muere en unos pocos centímetros.
3. Como $\sigma_t = 1.0\text{ cm}^{-1}$, los valores numéricos de $\nu$ coinciden con los de $\sigma_t\nu$.

---

## 5. Deducción de los Autovectores Angulares

Para cada autovalor se evalúa $v_m(\nu) = \dfrac{0.485}{1 - \mu_m/\nu}$ en las cuatro direcciones.

### 5.1 Modo fundamental ($\nu_1 = +3.3740096\text{ cm}$)

* **Para la dirección $m = 4$ ($\mu_4 = +0.8611363$):**
  $$v_4(\nu_1) = \frac{0.485}{1 - \frac{0.8611363}{3.3740096}} = \frac{0.485}{1 - 0.2552265} = \frac{0.485}{0.7447735} \approx 0.6512046$$

* **Para la dirección $m = 3$ ($\mu_3 = +0.3399810$):**
  $$v_3(\nu_1) = \frac{0.485}{1 - \frac{0.3399810}{3.3740096}} = \frac{0.485}{1 - 0.1007646} = \frac{0.485}{0.8992354} \approx 0.5393472$$

* **Para la dirección $m = 2$ ($\mu_2 = -0.3399810$):**
  $$v_2(\nu_1) = \frac{0.485}{1 + 0.1007646} = \frac{0.485}{1.1007646} \approx 0.4406028$$

* **Para la dirección $m = 1$ ($\mu_1 = -0.8611363$):**
  $$v_1(\nu_1) = \frac{0.485}{1 + 0.2552265} = \frac{0.485}{1.2552265} \approx 0.3863845$$

### 5.2 Modo transitorio ($\nu_2 = +0.5009792\text{ cm}$)

* **$m = 4$:** $\quad v_4(\nu_2) = \dfrac{0.485}{1 - \frac{0.8611363}{0.5009792}} = \dfrac{0.485}{1 - 1.7189217} = \dfrac{0.485}{-0.7189217} \approx -0.6746358$

* **$m = 3$:** $\quad v_3(\nu_2) = \dfrac{0.485}{1 - 0.6786316} = \dfrac{0.485}{0.3213684} \approx 1.5091784$

* **$m = 2$:** $\quad v_2(\nu_2) = \dfrac{0.485}{1 + 0.6786316} = \dfrac{0.485}{1.6786316} \approx 0.2889256$

* **$m = 1$:** $\quad v_1(\nu_2) = \dfrac{0.485}{1 + 1.7189217} = \dfrac{0.485}{2.7189217} \approx 0.1783805$

### 5.3 Tabla de autovectores

| $m$ | $\mu_m$ | $v_m(+\nu_1)$ | $v_m(-\nu_1)$ | $v_m(+\nu_2)$ | $v_m(-\nu_2)$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | $-0.8611363$ | $0.3863845$ | $0.6512046$ | $0.1783805$ | $-0.6746358$ |
| 2 | $-0.3399810$ | $0.4406028$ | $0.5393472$ | $0.2889256$ | $1.5091784$ |
| 3 | $+0.3399810$ | $0.5393472$ | $0.4406028$ | $1.5091784$ | $0.2889256$ |
| 4 | $+0.8611363$ | $0.6512046$ | $0.3863845$ | $-0.6746358$ | $0.1783805$ |

**Propiedad de simetría angular:** $v_m(-\nu_j) = v_{5-m}(\nu_j)$, es decir, el autovector del modo $-\nu_j$ es la imagen especular del modo $+\nu_j$ (la misma propiedad que en $S_2$).

**Verificación de la normalización:** para cada modo se cumple $\sum_m w_m v_m(\nu_j) = 1$:

$$\sum_m w_m v_m(\nu_1) = 0.3478548(0.3863845) + 0.6521452(0.4406028) + 0.6521452(0.5393472) + 0.3478548(0.6512046) = 1.0000$$
$$\sum_m w_m v_m(\nu_2) = 0.3478548(0.1783805) + 0.6521452(0.2889256) + 0.6521452(1.5091784) + 0.3478548(-0.6746358) = 1.0000$$

Nótese que el modo transitorio $\nu_2$ es fuertemente anisótropo e incluso **cambia de signo** entre las direcciones entrantes y la dirección saliente $m=4$. Esto es esperable: representa la corrección angular necesaria para ajustar el flujo impuesto en la frontera, imposible de reproducir con un único modo casi isótropo.

---

## 6. Superposición General y Aplicación de Condiciones de Borde

### 6.1 Solución general

La solución es la combinación lineal de los cuatro modos. Conviene escribirlos separando los que decaen desde la izquierda (autovalores $+\nu_j$) de los que decaen desde la derecha (autovalores $-\nu_j$):

$$\psi_m(x) = \sum_{j=1}^{2} A_j \, v_m(+\nu_j) \, e^{-x/\nu_j} + \sum_{j=1}^{2} B_j \, v_m(-\nu_j) \, e^{-(L-x)/\nu_j}$$

### 6.2 Los términos de la frontera derecha son despreciables

El medio tiene $100\text{ cm}$ de espesor, equivalente a $L/\nu_1 = 29.638$ longitudes de relajación fundamentales. Los factores de atenuación son:

$$e^{-L/\nu_1} = e^{-29.638} \approx 1.34 \times 10^{-13}, \qquad e^{-L/\nu_2} = e^{-199.61} \approx 2.05 \times 10^{-87}$$

Resolviendo el sistema completo de cuatro ecuaciones se obtiene $B_1 \approx -1.42\times10^{-13}$ y $B_2 \approx -1.27\times10^{-14}$, de modo que los modos provenientes de la derecha contribuyen con menos de $10^{-12}$ en todo el dominio, salvo en una capa de fracciones de centímetro junto a $x = 100\text{ cm}$. Por tanto, para imponer la condición de contorno izquierda pueden despreciarse, igual que en $S_2$ se tomó $B \approx 0$.

### 6.3 Sistema para las amplitudes $A_1$ y $A_2$

Aplicando la condición de flujo entrante unitario en $x = 0$ a las dos direcciones con $\mu_m > 0$ ($m = 3, 4$):

$$\psi_3(0) = A_1 v_3(\nu_1) + A_2 v_3(\nu_2) = 1.0$$
$$\psi_4(0) = A_1 v_4(\nu_1) + A_2 v_4(\nu_2) = 1.0$$

que con los valores de la tabla es el sistema $2\times2$:

$$0.5393472\, A_1 + 1.5091784\, A_2 = 1.0$$
$$0.6512046\, A_1 - 0.6746358\, A_2 = 1.0$$

Determinante:

$$\Delta = (0.5393472)(-0.6746358) - (1.5091784)(0.6512046) = -0.363857 - 0.982746 = -1.3466468$$

Regla de Cramer:

$$A_1 = \frac{(1.0)(-0.6746358) - (1.5091784)(1.0)}{\Delta} = \frac{-2.1838142}{-1.3466468} \approx 1.6216681$$

$$A_2 = \frac{(0.5393472)(1.0) - (0.6512046)(1.0)}{\Delta} = \frac{-0.1118574}{-1.3466468} \approx 0.0830637$$

Comparado con $S_2$ ($A = 1.7047318$, único), la amplitud del modo fundamental **disminuye** ($1.6216681$) porque parte del flujo incidente se describe ahora mediante el modo transitorio ($A_2 = 0.0830637$).

### 6.4 Verificación: el flujo escalar en la frontera es independiente de $N$

Un resultado notable se obtiene al sumar las dos amplitudes. Escribiendo $v_m(\nu_j) = \dfrac{(c/2)\nu_j}{\nu_j - \mu_m}$ y definiendo $S_j = A_j\nu_j$, la condición de contorno en $m = 3,4$ se escribe:

$$\sum_{j=1}^{2} \frac{S_j}{\nu_j - \mu_m} = \frac{2}{c}, \qquad m = 3, 4$$

Sea $R(x) = \sum_j \dfrac{S_j}{\nu_j - x}$. Entonces $R(\mu_3) = R(\mu_4) = 2/c$ y

$$R(x) - \frac{2}{c} = \frac{N(x)}{(\nu_1 - x)(\nu_2 - x)}$$

donde $N(x)$ es un polinomio de grado 2 con coeficiente principal $-2/c$ y raíces $\mu_3, \mu_4$:

$$N(x) = -\frac{2}{c}(x - \mu_3)(x - \mu_4)$$

Evaluando en $x = 0$ y usando $\mu_3\mu_4 = \mu_1\mu_2$:

$$R(0) = A_1 + A_2 = \frac{2}{c} \left[ 1 - \frac{\mu_1\mu_2}{\nu_1\nu_2} \right]$$

Como el producto de las raíces de $63\nu^4 - 733\nu^2 + 180 = 0$ es $\nu_1^2\nu_2^2 = 180/63 = 20/7$, resulta $\nu_1\nu_2 = \sqrt{20/7} = 1.6903085$, y con $\mu_1\mu_2 = 0.2927700$:

$$\frac{\mu_1\mu_2}{\nu_1\nu_2} = \sqrt{1 - c} = \sqrt{0.03} = 0.1732051$$

Por lo tanto:

$$\boxed{\;\phi(0) = A_1 + A_2 = \frac{2}{c}\Big(1 - \sqrt{1-c}\Big) = \frac{2}{1 + \sqrt{1-c}} = 1.7047318\;}$$

Este valor coincide **exactamente** con el obtenido en $S_2$ ($A = 1.7047318$): el flujo escalar en la superficie iluminada no depende del orden de la cuadratura. Lo que sí mejora con $N$ es la **distribución angular** del flujo saliente (y, con ella, la corriente y el albedo), no su valor integral.

