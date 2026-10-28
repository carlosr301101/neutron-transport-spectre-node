# Resolución Analítica del Problema de Transporte Neutrónico 1D ($S_2$) mediante Análisis Espectral

## 1. Definición y Datos del Problema

Se considera un medio homogéneo unidimensional (placa plana) sin fuentes internas fijas ($Q = 0$) ni fisión, con flujo neutrónico incidente por la frontera izquierda.

### Parámetros Físicos y Geométricos
* **Espesor de la placa ($L$):** $100\text{ cm}$ ($x \in [0, 100]$)
* **Sección eficaz total ($\sigma_t$):** $1.0\text{ cm}^{-1}$
* **Sección eficaz de dispersión isotrópica ($\sigma_{s0}$):** $0.97\text{ cm}^{-1}$
* **Razón de dispersión / Albedo de colisión ($c$):** $c = \frac{\sigma_{s0}}{\sigma_t} = \frac{0.97}{1.0} = 0.97$

### Condiciones de Contorno
* **Frontera Izquierda ($x = 0\text{ cm}$):** Flujo entrante conocido $\psi_m(0) = 1.0$ para todas las direcciones coincidentes con $\mu_m > 0$.
* **Frontera Derecha ($x = 100\text{ cm}$):** Condición de vacío $\psi_m(100) = 0.0$ para todas las direcciones $\mu_m < 0$.

---

## 2. Discretización Angular en Ordenadas Discretas ($S_2$)

Para la ordenación discreta $S_2$ ($N = 2$), los puntos de evaluación $\mu_m$ y los pesos $w_m$ provienen de la cuadratura de Gauss-Legendre (raíces del polinomio $P_2(\mu) = 0$):

* **Direcciones discretas:**
  $$\mu_1 = -\frac{1}{\sqrt{3}} \approx -0.5773503 \quad (\text{hacia la izquierda})$$
  $$\mu_2 = +\frac{1}{\sqrt{3}} \approx +0.5773503 \quad (\text{hacia la derecha})$$
* **Pesos de integración:** $w_1 = 1.0$, $w_2 = 1.0$ ($w_1 + w_2 = 2$).

El sistema acoplado de ecuaciones diferenciales ordinarias $S_2$ se expresa como:

$$\mu_1 \frac{d\psi_1(x)}{dx} + \sigma_t \psi_1(x) = \frac{\sigma_{s0}}{2} \Big[ \psi_1(x) + \psi_2(x) \Big]$$

$$\mu_2 \frac{d\psi_2(x)}{dx} + \sigma_t \psi_2(x) = \frac{\sigma_{s0}}{2} \Big[ \psi_1(x) + \psi_2(x) \Big]$$

donde el flujo escalar asociado a esta cuadratura es

$$\phi(x) = \sum_m w_m \psi_m(x) = 1.0 \cdot \psi_1(x) + 1.0 \cdot \psi_2(x) = \psi_1(x) + \psi_2(x)$$

ya que $w_1 = w_2 = 1.0$ ($\sum_m w_m = 2$, no $1$). Este factor es el que aparece en el término de fuente $\frac{\sigma_{s0}}{2}\phi(x)$.

---

## 3. Separación de Variables (Ansatz Espectral Modal)

Planteamos la solución mediante modales exponenciales homogéneos de la forma:

$$\psi_m(x) = v_m(\nu) e^{-x/\nu}, \quad m = 1, 2$$

donde:
* $\nu$ es el autovalor espectral (longitud de relajación o atenuación espacial).
* $v_m(\nu)$ representa la distribución angular del autovector correspondiente.

Sustituyendo el ansatz en el sistema de ecuaciones:

$$\left( \sigma_t - \frac{\mu_m}{\nu} \right) v_m(\nu) = \frac{\sigma_{s0}}{2} \underbrace{\big[ v_1(\nu) + v_2(\nu) \big]}_{c_0 = 1}$$

Normalizando el flujo escalar modal a la unidad ($v_1 + v_2 = 1$), despejamos el autovector angular $v_m(\nu)$:

$$v_m(\nu) = \frac{c/2}{1 - \frac{\mu_m}{\sigma_t \nu}}$$

---

## 4. Deducción Analítica de los Autovalores ($\nu$)

Para determinar los autovalores $\nu$, imponemos la condición de consistencia del flujo escalar modal $v_1(\nu) + v_2(\nu) = 1$:

$$\frac{c}{2} \left[ \frac{1}{1 - \frac{\mu_1}{\sigma_t \nu}} + \frac{1}{1 - \frac{\mu_2}{\sigma_t \nu}} \right] = 1$$

Dado que $\mu_1 = -\mu_2$ y $\sigma_t = 1.0\text{ cm}^{-1}$:

$$\frac{c}{2} \left[ \frac{1}{1 + \frac{\mu_2}{\nu}} + \frac{1}{1 - \frac{\mu_2}{\nu}} \right] = \frac{c}{2} \left[ \frac{2}{1 - \frac{\mu_2^2}{\nu^2}} \right] = \frac{c}{1 - \frac{\mu_2^2}{\nu^2}} = 1$$

Despejando $\nu^2$:

$$\nu^2 (1 - c) = \mu_2^2 \implies \nu^2 = \frac{\mu_2^2}{1 - c}$$

Sustituyendo los valores numéricos exactos ($\mu_2^2 = 1/3$ y $1 - c = 0.03 = 3/100$):

$$\nu^2 = \frac{1/3}{3/100} = \frac{100}{9} \implies \nu = \pm \frac{10}{3} \text{ cm}$$

Obtenemos dos autovalores simétricos:
$$\nu_+ = +\frac{10}{3}\text{ cm} \approx +3.3333\text{ cm}$$
$$\nu_- = -\frac{10}{3}\text{ cm} \approx -3.3333\text{ cm}$$

---

## 5. Deducción de los Autovectores Angulares

Sustituyendo $\nu_+ = +10/3\text{ cm}$ en la fórmula de $v_m(\nu)$:

* **Para la dirección derecha ($\mu_2 = +1/\sqrt{3}$):**
  $$v_{2,+} = \frac{0.485}{1 - \frac{1/\sqrt{3}}{10/3}} = \frac{0.485}{1 - \frac{\sqrt{3}}{10}} = \frac{4.85}{10 - \sqrt{3}} \approx 0.5866025$$

* **Para la dirección izquierda ($\mu_1 = -1/\sqrt{3}$):**
  $$v_{1,+} = \frac{0.485}{1 + \frac{1/\sqrt{3}}{10/3}} = \frac{0.485}{1 + \frac{\sqrt{3}}{10}} = \frac{4.85}{10 + \sqrt{3}} \approx 0.4133975$$

Por propiedades de simetría angular, para $\nu_- = -10/3\text{ cm}$:
$$v_{2,-} = v_{1,+} \approx 0.4133975, \quad v_{1,-} = v_{2,+} \approx 0.5866025$$

---

## 6. Superposición General y Aplicación de Condiciones de Borde

La solución general se expresa como combinación lineal de ambos modos:

$$\psi_1(x) = A \cdot v_{1,+} \, e^{-0.3 x} + B \cdot v_{1,-} \, e^{-0.3 (100 - x)}$$
$$\psi_2(x) = A \cdot v_{2,+} \, e^{-0.3 x} + B \cdot v_{2,-} \, e^{-0.3 (100 - x)}$$

Dado que el medio tiene un espesor de $100\text{ cm}$ (equivalente a $30$ longitudes de relajación $\nu$), la atenuación $e^{-30} \approx 9.36 \times 10^{-14} \approx 0$. La condición de contorno de vacío en la derecha ($\psi_1(100) = 0$) implica que $B \approx 0$.

Aplicando la condición de flujo incidente en la izquierda ($x = 0$):

$$\psi_2(0) = 1.0 \implies A \cdot v_{2,+} = 1.0 \implies A = \frac{1}{v_{2,+}} = \frac{10 - \sqrt{3}}{4.85} \approx 1.7047318$$

---

## 7. Resultados Finales y Análisis Físico

### Flujos Angulares Discretos
$$\psi_1(x) = 1.7047318 \times 0.4133975 \, e^{-0.3 x} \approx \mathbf{0.7047318 \, e^{-0.3 x}}$$
$$\psi_2(x) = 1.7047318 \times 0.5866025 \, e^{-0.3 x} = \mathbf{1.0000000 \, e^{-0.3 x}}$$

### Flujo Escalar Neutrónico $\phi(x)$
$$\phi(x) = \psi_1(x) + \psi_2(x) = 1.7047318 \Big( 0.4133975 + 0.5866025 \Big) e^{-0.3 x} = 1.7047318 \, e^{-0.3 x}$$
$$\mathbf{\phi(x) \approx 1.7047318 \, e^{-0.3 x}}$$

### Conclusiones Físicas
1. **Atenuación Espacial:** El flujo decrece de forma exponencial pura con un coeficiente constante de atenuación de $1/\nu = 0.3\text{ cm}^{-1}$.
2. **Respuesta de Reflejabilidad (Albedo):** En la frontera de entrada ($x = 0\text{ cm}$), el flujo emergente es $\psi_1(0) \approx 0.70473$ frente al flujo incidente $\psi_2(0) = 1.0$. El albedo se define como la razón de corrientes en la frontera, $J_{\text{out}}/J_{\text{in}} = \frac{w_1|\mu_1|\psi_1(0)}{w_2\mu_2\psi_2(0)} = \frac{\psi_1(0)}{\psi_2(0)} \approx 0.7047$, de modo que el **$70.47\%$ de los neutrones incidentes son reflejados** de vuelta hacia el exterior por el elevado valor de dispersión del medio ($c = 0.97$).
