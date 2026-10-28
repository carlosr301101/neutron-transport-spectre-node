# Resolução Analítica do Problema de Transporte de Nêutrons 1D ($S_2$) mediante Análise Espectral

## 1. Definição e Dados do Problema

Considera-se um meio homogêneo unidimensional (placa plana), sem fontes internas fixas ($Q = 0$) e sem fissão, com fluxo de nêutrons incidente pela fronteira esquerda.

### Parâmetros Físicos e Geométricos
* **Espessura da placa ($L$):** $100\text{ cm}$ ($x \in [0, 100]$)
* **Seção de choque total ($\sigma_t$):** $1.0\text{ cm}^{-1}$
* **Seção de choque de espalhamento isotrópico ($\sigma_{s0}$):** $0.97\text{ cm}^{-1}$
* **Razão de espalhamento / Albedo de colisão ($c$):** $c = \frac{\sigma_{s0}}{\sigma_t} = \frac{0.97}{1.0} = 0.97$

### Condições de Contorno
* **Fronteira Esquerda ($x = 0\text{ cm}$):** Fluxo entrante conhecido $\psi_m(0) = 1.0$ para todas as direções coincidentes com $\mu_m > 0$.
* **Fronteira Direita ($x = 100\text{ cm}$):** Condição de vácuo $\psi_m(100) = 0.0$ para todas as direções $\mu_m < 0$.

---

## 2. Discretização Angular em Ordenadas Discretas ($S_2$)

Para a ordenação discreta $S_2$ ($N = 2$), os pontos de avaliação $\mu_m$ e os pesos $w_m$ provêm da quadratura de Gauss-Legendre (raízes do polinômio $P_2(\mu) = 0$):

* **Direções discretas:**
  $$\mu_1 = -\frac{1}{\sqrt{3}} \approx -0.5773503 \quad (\text{para a esquerda})$$
  $$\mu_2 = +\frac{1}{\sqrt{3}} \approx +0.5773503 \quad (\text{para a direita})$$
* **Pesos de integração:** $w_1 = 1.0$, $w_2 = 1.0$ ($w_1 + w_2 = 2$).

O sistema acoplado de equações diferenciais ordinárias $S_2$ é expresso como:

$$\mu_1 \frac{d\psi_1(x)}{dx} + \sigma_t \psi_1(x) = \frac{\sigma_{s0}}{2} \Big[ \psi_1(x) + \psi_2(x) \Big]$$

$$\mu_2 \frac{d\psi_2(x)}{dx} + \sigma_t \psi_2(x) = \frac{\sigma_{s0}}{2} \Big[ \psi_1(x) + \psi_2(x) \Big]$$

em que o fluxo escalar associado a esta quadratura é

$$\phi(x) = \sum_m w_m \psi_m(x) = 1.0 \cdot \psi_1(x) + 1.0 \cdot \psi_2(x) = \psi_1(x) + \psi_2(x)$$

visto que $w_1 = w_2 = 1.0$ ($\sum_m w_m = 2$, e não $1$). Esse fator é o que aparece no termo de fonte $\frac{\sigma_{s0}}{2}\phi(x)$.

---

## 3. Separação de Variáveis (Ansatz Espectral Modal)

Propomos a solução por meio de modos exponenciais homogêneos da forma:

$$\psi_m(x) = v_m(\nu) e^{-x/\nu}, \quad m = 1, 2$$

em que:
* $\nu$ é o autovalor espectral (comprimento de relaxação ou atenuação espacial).
* $v_m(\nu)$ representa a distribuição angular do autovetor correspondente.

Substituindo o ansatz no sistema de equações:

$$\left( \sigma_t - \frac{\mu_m}{\nu} \right) v_m(\nu) = \frac{\sigma_{s0}}{2} \underbrace{\big[ v_1(\nu) + v_2(\nu) \big]}_{c_0 = 1}$$

Normalizando o fluxo escalar modal à unidade ($v_1 + v_2 = 1$), isolamos o autovetor angular $v_m(\nu)$:

$$v_m(\nu) = \frac{c/2}{1 - \frac{\mu_m}{\sigma_t \nu}}$$

---

## 4. Dedução Analítica dos Autovalores ($\nu$)

Para determinar os autovalores $\nu$, impomos a condição de consistência do fluxo escalar modal $v_1(\nu) + v_2(\nu) = 1$:

$$\frac{c}{2} \left[ \frac{1}{1 - \frac{\mu_1}{\sigma_t \nu}} + \frac{1}{1 - \frac{\mu_2}{\sigma_t \nu}} \right] = 1$$

Visto que $\mu_1 = -\mu_2$ e $\sigma_t = 1.0\text{ cm}^{-1}$:

$$\frac{c}{2} \left[ \frac{1}{1 + \frac{\mu_2}{\nu}} + \frac{1}{1 - \frac{\mu_2}{\nu}} \right] = \frac{c}{2} \left[ \frac{2}{1 - \frac{\mu_2^2}{\nu^2}} \right] = \frac{c}{1 - \frac{\mu_2^2}{\nu^2}} = 1$$

Isolando $\nu^2$:

$$\nu^2 (1 - c) = \mu_2^2 \implies \nu^2 = \frac{\mu_2^2}{1 - c}$$

Substituindo os valores numéricos exatos ($\mu_2^2 = 1/3$ e $1 - c = 0.03 = 3/100$):

$$\nu^2 = \frac{1/3}{3/100} = \frac{100}{9} \implies \nu = \pm \frac{10}{3} \text{ cm}$$

Obtemos dois autovalores simétricos:
$$\nu_+ = +\frac{10}{3}\text{ cm} \approx +3.3333\text{ cm}$$
$$\nu_- = -\frac{10}{3}\text{ cm} \approx -3.3333\text{ cm}$$

---

## 5. Dedução dos Autovetores Angulares

Substituindo $\nu_+ = +10/3\text{ cm}$ na fórmula de $v_m(\nu)$:

* **Para a direção direita ($\mu_2 = +1/\sqrt{3}$):**
  $$v_{2,+} = \frac{0.485}{1 - \frac{1/\sqrt{3}}{10/3}} = \frac{0.485}{1 - \frac{\sqrt{3}}{10}} = \frac{4.85}{10 - \sqrt{3}} \approx 0.5866025$$

* **Para a direção esquerda ($\mu_1 = -1/\sqrt{3}$):**
  $$v_{1,+} = \frac{0.485}{1 + \frac{1/\sqrt{3}}{10/3}} = \frac{0.485}{1 + \frac{\sqrt{3}}{10}} = \frac{4.85}{10 + \sqrt{3}} \approx 0.4133975$$

Por propriedades de simetria angular, para $\nu_- = -10/3\text{ cm}$:
$$v_{2,-} = v_{1,+} \approx 0.4133975, \quad v_{1,-} = v_{2,+} \approx 0.5866025$$

---

## 6. Superposição Geral e Aplicação das Condições de Contorno

A solução geral é expressa como combinação linear dos dois modos:

$$\psi_1(x) = A \cdot v_{1,+} \, e^{-0.3 x} + B \cdot v_{1,-} \, e^{-0.3 (100 - x)}$$
$$\psi_2(x) = A \cdot v_{2,+} \, e^{-0.3 x} + B \cdot v_{2,-} \, e^{-0.3 (100 - x)}$$

Visto que o meio tem uma espessura de $100\text{ cm}$ (equivalente a $30$ comprimentos de relaxação $\nu$), a atenuação $e^{-30} \approx 9.36 \times 10^{-14} \approx 0$. A condição de contorno de vácuo à direita ($\psi_1(100) = 0$) acarreta $B \approx 0$.

Aplicando a condição de fluxo incidente à esquerda ($x = 0$):

$$\psi_2(0) = 1.0 \implies A \cdot v_{2,+} = 1.0 \implies A = \frac{1}{v_{2,+}} = \frac{10 - \sqrt{3}}{4.85} \approx 1.7047318$$

---

## 7. Resultados Finais

### Fluxos Angulares Discretos
$$\psi_1(x) = 1.7047318 \times 0.4133975 \, e^{-0.3 x} \approx \mathbf{0.7047318 \, e^{-0.3 x}}$$
$$\psi_2(x) = 1.7047318 \times 0.5866025 \, e^{-0.3 x} = \mathbf{1.0000000 \, e^{-0.3 x}}$$

### Fluxo Escalar de Nêutrons $\phi(x)$
$$\phi(x) = \psi_1(x) + \psi_2(x) = 1.7047318 \Big( 0.4133975 + 0.5866025 \Big) e^{-0.3 x} = 1.7047318 \, e^{-0.3 x}$$
$$\mathbf{\phi(x) \approx 1.7047318 \, e^{-0.3 x}}$$

