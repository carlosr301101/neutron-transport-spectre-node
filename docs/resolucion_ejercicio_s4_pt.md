# Resolução Analítica do Problema de Transporte de Nêutrons 1D ($S_4$) mediante Análise Espectral

## 1. Definição e Dados do Problema

Considera-se o mesmo meio homogêneo unidimensional (placa plana) do exercício $S_2$: sem fontes internas fixas ($Q = 0$), sem fissão, com fluxo de nêutrons incidente pela fronteira esquerda. A única diferença é o **aumento da ordem da quadratura angular**, de $S_2$ ($N = 2$) para $S_4$ ($N = 4$), o que duplica o número de direções e, portanto, o número de modos espectrais.

### Parâmetros Físicos e Geométricos
* **Espessura da placa ($L$):** $100\text{ cm}$ ($x \in [0, 100]$)
* **Seção de choque total ($\sigma_t$):** $1.0\text{ cm}^{-1}$
* **Seção de choque de espalhamento isotrópico ($\sigma_{s0}$):** $0.97\text{ cm}^{-1}$
* **Razão de espalhamento / Albedo de colisão ($c$):** $c = \frac{\sigma_{s0}}{\sigma_t} = \frac{0.97}{1.0} = 0.97$

### Condições de Contorno
* **Fronteira Esquerda ($x = 0\text{ cm}$):** Fluxo entrante conhecido $\psi_m(0) = 1.0$ para todas as direções com $\mu_m > 0$ (ou seja, $m = 3, 4$).
* **Fronteira Direita ($x = 100\text{ cm}$):** Condição de vácuo $\psi_m(100) = 0.0$ para todas as direções com $\mu_m < 0$ (ou seja, $m = 1, 2$).

Em $S_4$ há **duas** direções entrantes e **duas** emergentes em cada fronteira, de modo que se impõem quatro condições de contorno para determinar os quatro modos da solução.

---

## 2. Discretização Angular em Ordenadas Discretas ($S_4$)

Para a ordenação discreta $S_4$ ($N = 4$), os pontos de avaliação $\mu_m$ e os pesos $w_m$ provêm da quadratura de Gauss-Legendre (raízes do polinômio $P_4(\mu) = 0$), cujos valores são tomados do arquivo `cuadraturas.py`:

$$\{\mu_k\}_{\mu>0} = \{\,0.861136311594053,\; 0.339981043584856\,\}$$
$$\{w_k\}_{\mu>0} = \{\,0.347854845137454,\; 0.652145154862546\,\}$$

Na forma fechada exata, são:

$$\mu = \pm \sqrt{\frac{3 \mp 2\sqrt{6/5}}{7}}, \qquad w = \frac{18 \pm \sqrt{30}}{36}$$

Como em $S_2$, o conjunto de direções é simétrico ($\mu_{m'} = -\mu_m$ para o índice espelho) e os pesos são iguais para cada par $\pm\mu_k$. Reordenando os índices da esquerda para a direita, tem-se:

| $m$ | $\mu_m$ | $w_m$ | Sentido |
|:---:|:---:|:---:|:---|
| 1 | $-0.8611363$ | $0.3478548$ | para a esquerda |
| 2 | $-0.3399810$ | $0.6521452$ | para a esquerda |
| 3 | $+0.3399810$ | $0.6521452$ | para a direita |
| 4 | $+0.8611363$ | $0.3478548$ | para a direita |

* **Pesos de integração:** $w_1 + w_2 + w_3 + w_4 = 2.0$ (igual a $S_2$).

O sistema acoplado de quatro equações diferenciais ordinárias $S_4$ é:

$$\mu_m \frac{d\psi_m(x)}{dx} + \sigma_t \psi_m(x) = \frac{\sigma_{s0}}{2} \sum_{n=1}^{4} w_n \psi_n(x), \quad m = 1, 2, 3, 4$$

ou seja, explicitamente:

$$\mu_1 \frac{d\psi_1}{dx} + \sigma_t \psi_1 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_2 \frac{d\psi_2}{dx} + \sigma_t \psi_2 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_3 \frac{d\psi_3}{dx} + \sigma_t \psi_3 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$
$$\mu_4 \frac{d\psi_4}{dx} + \sigma_t \psi_4 = \frac{\sigma_{s0}}{2} \Big[ w_1\psi_1 + w_2\psi_2 + w_3\psi_3 + w_4\psi_4 \Big]$$

em que o fluxo escalar associado a esta quadratura é

$$\phi(x) = \sum_{m=1}^{4} w_m \psi_m(x), \qquad \sum_{m=1}^{4} w_m = 2$$

Esse é o fator que aparece no termo de fonte $\frac{\sigma_{s0}}{2}\phi(x)$; como em $S_2$, não há um fator $1/2$ adicional.

---

## 3. Separação de Variáveis (Ansatz Espectral Modal)

Propomos a solução por meio de modos exponenciais homogêneos da forma:

$$\psi_m(x) = v_m(\nu) e^{-x/\nu}, \quad m = 1, 2, 3, 4$$

em que:
* $\nu$ é o autovalor espectral (comprimento de relaxação ou atenuação espacial).
* $v_m(\nu)$ é a componente $m$-ésima do autovetor angular correspondente.

Substituindo o ansatz no sistema de equações e usando que $\frac{d\psi_m}{dx} = -\frac{1}{\nu} v_m e^{-x/\nu}$:

$$\left( \sigma_t - \frac{\mu_m}{\nu} \right) v_m(\nu) = \frac{\sigma_{s0}}{2} \underbrace{\sum_{n=1}^{4} w_n v_n(\nu)}_{\text{fluxo escalar modal}}$$

Normalizando o fluxo escalar modal à unidade ($\sum_n w_n v_n(\nu) = 1$) e lembrando que $\sigma_{s0}/\sigma_t = c$, isola-se o autovetor angular:

$$v_m(\nu) = \frac{c/2}{1 - \frac{\mu_m}{\sigma_t \nu}} = \frac{0.485}{1 - \frac{\mu_m}{\nu}} \qquad (\sigma_t = 1.0\text{ cm}^{-1})$$

Diferentemente de $S_2$ (um único par de autovalores $\pm\nu$), em $S_4$ a condição de consistência leva a uma equação de **quarto grau**, de modo que existirão **dois pares** de autovalores $\pm\nu_1$ e $\pm\nu_2$ e, com eles, **quatro modos** espectrais.

---

## 4. Dedução Analítica dos Autovalores ($\nu$)

### 4.1 Condição de consistência

Impondo $\sum_n w_n v_n(\nu) = 1$ na fórmula do autovetor:

$$\frac{c}{2} \sum_{m=1}^{4} \frac{w_m}{1 - \frac{\mu_m}{\sigma_t \nu}} = 1 \quad \Longrightarrow \quad \sum_{m=1}^{4} \frac{w_m}{1 - \frac{\mu_m}{\sigma_t \nu}} = \frac{2}{c}$$

### 4.2 Redução por simetria

Como $\mu_{1,2} = -\mu_{3,4}$ com pesos idênticos, agrupam-se os pares $\pm\mu_k$ ($k=1,2$):

$$\sum_{k=1}^{2} w_k \left[ \frac{1}{1 - \frac{\mu_k}{\sigma_t \nu}} + \frac{1}{1 + \frac{\mu_k}{\sigma_t \nu}} \right] = \frac{2}{c}$$

A soma entre colchetes é $\dfrac{2}{1 - \left(\frac{\mu_k}{\sigma_t \nu}\right)^2}$, de modo que

$$\sum_{k=1}^{2} \frac{w_k}{1 - \frac{\mu_k^2}{(\sigma_t \nu)^2}} = \frac{1}{c}$$

e com $\sigma_t = 1.0\text{ cm}^{-1}$ (numericamente $\sigma_t\nu = \nu$):

$$\frac{w_1 \nu^2}{\nu^2 - \mu_1^2} + \frac{w_2 \nu^2}{\nu^2 - \mu_2^2} = \frac{1}{c}$$

### 4.3 Equação característica

Multiplicando por $(\nu^2 - \mu_1^2)(\nu^2 - \mu_2^2)$ e usando $w_1 + w_2 = 1$:

$$w_1 \nu^2 (\nu^2 - \mu_2^2) + w_2 \nu^2 (\nu^2 - \mu_1^2) = \frac{1}{c} (\nu^2 - \mu_1^2)(\nu^2 - \mu_2^2)$$

$$\nu^4 - Q \nu^2 = \frac{1}{c} \Big( \nu^4 - P \nu^2 + R \Big)$$

em que se definiram as constantes (que dependem apenas da quadratura):

$$P = \mu_1^2 + \mu_2^2, \qquad Q = w_1\mu_2^2 + w_2\mu_1^2, \qquad R = \mu_1^2\mu_2^2$$

Reagrupando, obtém-se uma **equação biquadrada** em $\nu$:

$$(1 - c)\nu^4 + (cQ - P)\nu^2 + R = 0$$

### 4.4 Valores numéricos exatos

Com $\mu_1^2 = 0.7415557$, $\mu_2^2 = 0.1155871$, $w_1 = 0.3478548$, $w_2 = 0.6521452$:

$$P = \mu_1^2 + \mu_2^2 = 0.8571429 = \frac{6}{7}$$

$$Q = w_1\mu_2^2 + w_2\mu_1^2 = 0.5238095 = \frac{11}{21}$$

$$R = \mu_1^2\mu_2^2 = 0.0857143 = \frac{3}{35}$$

$$cQ - P = 0.97 \cdot \frac{11}{21} - \frac{6}{7} = \frac{10.67}{21} - \frac{18}{21} = -\frac{733}{2100} \approx -0.3490476$$

Substituindo na equação característica:

$$\frac{3}{100}\nu^4 - \frac{733}{2100}\nu^2 + \frac{3}{35} = 0$$

Multiplicando tudo por $2100$ obtém-se a **equação característica exata de $S_4$**:

$$\boxed{\;63\nu^4 - 733\nu^2 + 180 = 0\;}$$

(Compare-se com $S_2$, em que a condição se reduzia à simples $\nu^2 = 100/9$.)

### 4.5 Raízes

Com $y = \nu^2$, a equação se reduz à quadrática $63y^2 - 733y + 180 = 0$:

$$y = \frac{733 \pm \sqrt{733^2 - 4 \cdot 63 \cdot 180}}{2 \cdot 63} = \frac{733 \pm \sqrt{491929}}{126}$$

$$\sqrt{491929} = 701.3765037$$

$$y_1 = \nu_1^2 = \frac{733 + 701.3765037}{126} = 11.3839405 \quad \Longrightarrow \quad \nu_1 = 3.3740096\text{ cm}$$
$$y_2 = \nu_2^2 = \frac{733 - 701.3765037}{126} = 0.2509801 \quad \Longrightarrow \quad \nu_2 = 0.5009792\text{ cm}$$

Obtemos assim **quatro** autovalores simétricos (dois pares):

| Autovalor | $\nu_j$ | $1/\nu_j$ ($\text{cm}^{-1}$) | Tipo de modo |
|:---:|:---:|:---:|:---|
| $\nu_1$ | $+3.3740096$ | $0.2963833$ | fundamental (decai para a direita) |
| $-\nu_1$ | $-3.3740096$ | $-0.2963833$ | fundamental (cresce para a direita) |
| $\nu_2$ | $+0.5009792$ | $1.9960910$ | transitório (camada limite esquerda) |
| $-\nu_2$ | $-0.5009792$ | $-1.9960910$ | transitório (camada limite direita) |

**Observações:**
1. O autovalor fundamental $\nu_1 = 3.3740096\text{ cm}$ é maior que o de $S_2$ ($10/3 = 3.3333\text{ cm}$): ao refinar a quadratura, o comprimento de relaxação aumenta, ou seja, o fluxo decai mais lentamente.
2. O segundo autovalor $\nu_2 = 0.5009792\text{ cm}$ é muito menor que a espessura: define um modo transitório que morre em poucos centímetros.
3. Como $\sigma_t = 1.0\text{ cm}^{-1}$, os valores numéricos de $\nu$ coincidem com os de $\sigma_t\nu$.

---

## 5. Dedução dos Autovetores Angulares

Para cada autovalor avalia-se $v_m(\nu) = \dfrac{0.485}{1 - \mu_m/\nu}$ nas quatro direções.

### 5.1 Modo fundamental ($\nu_1 = +3.3740096\text{ cm}$)

* **Para a direção $m = 4$ ($\mu_4 = +0.8611363$):**
  $$v_4(\nu_1) = \frac{0.485}{1 - \frac{0.8611363}{3.3740096}} = \frac{0.485}{1 - 0.2552265} = \frac{0.485}{0.7447735} \approx 0.6512046$$

* **Para a direção $m = 3$ ($\mu_3 = +0.3399810$):**
  $$v_3(\nu_1) = \frac{0.485}{1 - \frac{0.3399810}{3.3740096}} = \frac{0.485}{1 - 0.1007646} = \frac{0.485}{0.8992354} \approx 0.5393472$$

* **Para a direção $m = 2$ ($\mu_2 = -0.3399810$):**
  $$v_2(\nu_1) = \frac{0.485}{1 + 0.1007646} = \frac{0.485}{1.1007646} \approx 0.4406028$$

* **Para a direção $m = 1$ ($\mu_1 = -0.8611363$):**
  $$v_1(\nu_1) = \frac{0.485}{1 + 0.2552265} = \frac{0.485}{1.2552265} \approx 0.3863845$$

### 5.2 Modo transitório ($\nu_2 = +0.5009792\text{ cm}$)

* **$m = 4$:** $\quad v_4(\nu_2) = \dfrac{0.485}{1 - \frac{0.8611363}{0.5009792}} = \dfrac{0.485}{1 - 1.7189217} = \dfrac{0.485}{-0.7189217} \approx -0.6746358$

* **$m = 3$:** $\quad v_3(\nu_2) = \dfrac{0.485}{1 - 0.6786316} = \dfrac{0.485}{0.3213684} \approx 1.5091784$

* **$m = 2$:** $\quad v_2(\nu_2) = \dfrac{0.485}{1 + 0.6786316} = \dfrac{0.485}{1.6786316} \approx 0.2889256$

* **$m = 1$:** $\quad v_1(\nu_2) = \dfrac{0.485}{1 + 1.7189217} = \dfrac{0.485}{2.7189217} \approx 0.1783805$

### 5.3 Tabela de autovetores

| $m$ | $\mu_m$ | $v_m(+\nu_1)$ | $v_m(-\nu_1)$ | $v_m(+\nu_2)$ | $v_m(-\nu_2)$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | $-0.8611363$ | $0.3863845$ | $0.6512046$ | $0.1783805$ | $-0.6746358$ |
| 2 | $-0.3399810$ | $0.4406028$ | $0.5393472$ | $0.2889256$ | $1.5091784$ |
| 3 | $+0.3399810$ | $0.5393472$ | $0.4406028$ | $1.5091784$ | $0.2889256$ |
| 4 | $+0.8611363$ | $0.6512046$ | $0.3863845$ | $-0.6746358$ | $0.1783805$ |

**Propriedade de simetria angular:** $v_m(-\nu_j) = v_{5-m}(\nu_j)$, ou seja, o autovetor do modo $-\nu_j$ é a imagem especular do modo $+\nu_j$ (a mesma propriedade que em $S_2$).

**Verificação da normalização:** para cada modo vale $\sum_m w_m v_m(\nu_j) = 1$:

$$\sum_m w_m v_m(\nu_1) = 0.3478548(0.3863845) + 0.6521452(0.4406028) + 0.6521452(0.5393472) + 0.3478548(0.6512046) = 1.0000$$
$$\sum_m w_m v_m(\nu_2) = 0.3478548(0.1783805) + 0.6521452(0.2889256) + 0.6521452(1.5091784) + 0.3478548(-0.6746358) = 1.0000$$

Note-se que o modo transitório $\nu_2$ é fortemente anisotrópico e até **muda de sinal** entre as direções entrantes e a direção emergente $m=4$. Isso é esperado: representa a correção angular necessária para ajustar o fluxo imposto na fronteira, impossível de reproduzir com um único modo quase isotrópico.

---

## 6. Superposição Geral e Aplicação das Condições de Contorno

### 6.1 Solução geral

A solução é a combinação linear dos quatro modos. Convém escrevê-los separando os que decaem a partir da esquerda (autovalores $+\nu_j$) dos que decaem a partir da direita (autovalores $-\nu_j$):

$$\psi_m(x) = \sum_{j=1}^{2} A_j \, v_m(+\nu_j) \, e^{-x/\nu_j} + \sum_{j=1}^{2} B_j \, v_m(-\nu_j) \, e^{-(L-x)/\nu_j}$$

### 6.2 Os termos da fronteira direita são desprezíveis

O meio tem $100\text{ cm}$ de espessura, equivalente a $L/\nu_1 = 29.638$ comprimentos de relaxação fundamentais. Os fatores de atenuação são:

$$e^{-L/\nu_1} = e^{-29.638} \approx 1.34 \times 10^{-13}, \qquad e^{-L/\nu_2} = e^{-199.61} \approx 2.05 \times 10^{-87}$$

Resolvendo o sistema completo de quatro equações obtêm-se $B_1 \approx -1.42\times10^{-13}$ e $B_2 \approx -1.27\times10^{-14}$, de modo que os modos provenientes da direita contribuem com menos de $10^{-12}$ em todo o domínio, exceto numa camada de frações de centímetro junto a $x = 100\text{ cm}$. Portanto, para impor a condição de contorno esquerda eles podem ser desprezados, assim como em $S_2$ tomou-se $B \approx 0$.

### 6.3 Sistema para as amplitudes $A_1$ e $A_2$

Aplicando a condição de fluxo entrante unitário em $x = 0$ às duas direções com $\mu_m > 0$ ($m = 3, 4$):

$$\psi_3(0) = A_1 v_3(\nu_1) + A_2 v_3(\nu_2) = 1.0$$
$$\psi_4(0) = A_1 v_4(\nu_1) + A_2 v_4(\nu_2) = 1.0$$

que, com os valores da tabela, é o sistema $2\times2$:

$$0.5393472\, A_1 + 1.5091784\, A_2 = 1.0$$
$$0.6512046\, A_1 - 0.6746358\, A_2 = 1.0$$

Determinante:

$$\Delta = (0.5393472)(-0.6746358) - (1.5091784)(0.6512046) = -0.363857 - 0.982746 = -1.3466468$$

Regra de Cramer:

$$A_1 = \frac{(1.0)(-0.6746358) - (1.5091784)(1.0)}{\Delta} = \frac{-2.1838142}{-1.3466468} \approx 1.6216681$$

$$A_2 = \frac{(0.5393472)(1.0) - (0.6512046)(1.0)}{\Delta} = \frac{-0.1118574}{-1.3466468} \approx 0.0830637$$

Comparado com $S_2$ ($A = 1.7047318$, único), a amplitude do modo fundamental **diminui** ($1.6216681$) porque parte do fluxo incidente é agora descrita pelo modo transitório ($A_2 = 0.0830637$).

### 6.4 Verificação: o fluxo escalar na fronteira é independente de $N$

Um resultado notável é obtido ao somar as duas amplitudes. Escrevendo $v_m(\nu_j) = \dfrac{(c/2)\nu_j}{\nu_j - \mu_m}$ e definindo $S_j = A_j\nu_j$, a condição de contorno em $m = 3,4$ escreve-se:

$$\sum_{j=1}^{2} \frac{S_j}{\nu_j - \mu_m} = \frac{2}{c}, \qquad m = 3, 4$$

Seja $R(x) = \sum_j \dfrac{S_j}{\nu_j - x}$. Então $R(\mu_3) = R(\mu_4) = 2/c$ e

$$R(x) - \frac{2}{c} = \frac{N(x)}{(\nu_1 - x)(\nu_2 - x)}$$

em que $N(x)$ é um polinômio de grau 2 com coeficiente principal $-2/c$ e raízes $\mu_3, \mu_4$:

$$N(x) = -\frac{2}{c}(x - \mu_3)(x - \mu_4)$$

Avaliando em $x = 0$ e usando $\mu_3\mu_4 = \mu_1\mu_2$:

$$R(0) = A_1 + A_2 = \frac{2}{c} \left[ 1 - \frac{\mu_1\mu_2}{\nu_1\nu_2} \right]$$

Como o produto das raízes de $63\nu^4 - 733\nu^2 + 180 = 0$ é $\nu_1^2\nu_2^2 = 180/63 = 20/7$, resulta $\nu_1\nu_2 = \sqrt{20/7} = 1.6903085$ e, com $\mu_1\mu_2 = 0.2927700$:

$$\frac{\mu_1\mu_2}{\nu_1\nu_2} = \sqrt{1 - c} = \sqrt{0.03} = 0.1732051$$

Portanto:

$$\boxed{\;\phi(0) = A_1 + A_2 = \frac{2}{c}\Big(1 - \sqrt{1-c}\Big) = \frac{2}{1 + \sqrt{1-c}} = 1.7047318\;}$$

Esse valor coincide **exatamente** com o obtido em $S_2$ ($A = 1.7047318$): o fluxo escalar na superfície iluminada não depende da ordem da quadratura. O que de fato melhora com $N$ é a **distribuição angular** do fluxo emergente (e, com ela, a corrente e o albedo), não seu valor integral.

