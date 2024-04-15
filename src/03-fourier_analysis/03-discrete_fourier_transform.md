## Diskrete Fourier-Transformation

Wenn wir die Fourier-Transformation auf ein Signal am Computer durchführen
wollen, gibt es zwei Probleme in 
Gl. {{eqref: eq:fourier_transform_backward}}
(bzw. auch in Gleichungen anderer Normierungen).
Einerseits ist die Funktion $f(t)$ an (überabzählbar) unendlich vielen 
Stellen definiert und andererseits muss das Integral von 
$-\infty$ bis $\infty$ berechnet werden. Der Computer kann aber nur endlich
viele Werte speichern und berechnen. Deshalb müssen wir die Funktion
diskretisieren und nur über ein beschränktes Intervall integrieren.

### Theoretische Grundlagen
Wir wählen ein gleichmäßiges Zeitgrid 
$\{t_n\}_ {n = 0,\cdots N-1} = \{0, 1, \cdots, N-1\}$ mit Abstand 
$\Delta t = 1$, also $t_n = n$. Möchte man die Funktion an einem anderen
gleichmäßigem Grid evaluieren, kann das durch Verschiebung und Skalierung
der Funktion realisiert werden. Damit wird 
Gl. {{eqref: eq:fourier_transform_backward}} zu
$$
  F^{\mathrm{b}}_ k = \sum_{n=0}^{N-1} f(t_n) \eu^{-\iu \omega_k t_n} \Delta t
  = \sum_{n=0}^{N-1} f_n \eu^{-\iu \omega_k n}\,,
$$
mit $f_n = f(t_n)$, wobei wir die Ersetzung 
$\int \cdot\ \du t \to \sum \cdot\ \Delta t$ vorgenommen haben.

Weil wir die Funktion nur im begrenzten Intervall $[0, N-1]$ kennen, ist die
niedrigste Frequenz abgesehen von $0$ die Grundfrequenz $\omega_1 = 2\pi/N$.
Damit gilt $\omega_k = 2\pi k/N$ für $k = 0, \cdots, N-1$.
Die *diskrete Fourier-Transformation mit Rückwärtsnormierung* ist also durch
$$
  F^{\mathrm{b}}_ k = \sum_{n=0}^{N-1} f_n \eu^{-\iu 2\pi \frac{k}{N} n}
  {{numeq}}{eq:dft_backward}
$$
gegeben. 

Man kann sich leicht überzeugen, dass bei der Rücktransformation in diesem
Fall einen Normierungsfaktor von $1/N$ notwendig ist, also
$$
  f_n = \frac{1}{N} \sum_{k=0}^{N-1} F^{\mathrm{b}}_ k \eu^{\iu 2\pi \frac{k}{N} n}\,.
  {{numeq}}{eq:idft_backward}
$$
Gl. {{eqref: eq:idft_backward}} beschreibt die 
*inverse diskrete Fourier-Transformation*. Genau so wie bei der 
kontinuierlichen Fourier-Transformation kann man hier auch 
Vorwärtsnormierung und symmetrische Normierung definieren. Die genauen 
Formeln sind dann trivial und werden hier nicht explizit aufgeführt.
Da wir in diesem Kurs haupsächlich mit der Rückwärtsnormierung arbeiten
werden, einigen wir uns darauf, dass $F_k = F^{\mathrm{b}}_ k$ ist, wenn
nicht explizit anders angegeben.

Obwohl auf dem ersten Blick wir Kreisfrequenzen bis zu $2\pi \frac{N-1}{N}$
analysieren können, bringt die Diskretisierung des ursprünglichen Signals
eine Schwierigkeit mit sich. Um dieses Problem zu verstehen, betrachten wir
das Signal $f(t) = \sin(2\pi \frac{48}{50} t)$, welches eine Kreisfrequenz
$2\pi \frac{48}{50}$ besitzt und korrespondiert zu $k = 48$ bei $N = 50$:
![DFT Sampling mit Kreisfrequenz 0.96 * 2pi](../assets/figures/03-fourier_analysis/dft_sampling_figures/0048.png)
Obwohl das stetige Signal (blau) eine höhere Frequenz hat, kann diese durch 
die grobe Diskretisierung (orange) nicht dargestellt werden. 
Das diskretisierte Signal suggeriert eine viel niedrigere Frequenz, welche
in diesem Fall $2\pi \frac{2}{50}$ ist. 

Tatsächlich können die Kreisfrequenzen über $\pi$ nicht mehr erhalten werden.
Die folgende Animation zeigt das diskretizierte Signal für $k = 0, \cdots, 49$:
<p>
  <img src="../assets/figures/03-fourier_analysis/dft_sampling_figures/dft_sampling.gif" alt="DFT Sampling animation" />
</p>
Es kann gezeigt werden, dass die Korrespondenz
$$
  \omega = \begin{cases}
    \frac{2\pi}{N\cdot \Delta t} 
      (0, 1, \cdots, \frac{N}{2} - 1, -\frac{N}{2}, \cdots, -1)^\intercal & \text{für } N \text{ gerade} \\
    \frac{2\pi}{N\cdot \Delta t} 
      (0, 1, \cdots, \frac{N-1}{2}, -\frac{N-1}{2}, \cdots, -1)^\intercal & \text{für } N \text{ ungerade}
  \end{cases}
  {{numeq}}{eq:dft_frequencies}
$$

gilt. Das bedeutet, dass die obere Hälfte des Spektrums nicht die höheren
Kreisfrequenzen repräsentiert, sondern die negativen Frequenzen. 
Der Faktor $1/\Delta t$ in Gl. {{eqref: eq:dft_frequencies}} 
berücksichtigt die Skalierung des Zeitgrids im Fall $\Delta t \neq 1$.

Mit den Gleichungen {{eqref: eq:dft_backward}} und 
{{eqref: eq:dft_frequencies}} lässt sich die 
diskrete Fourier-Transformation implementieren.

### Implementierung
Wir implementieren die diskrete Fourier-Transformation am Beispiel des
Fourier-Transformations-Infrarot-Spektroskopie (FTIR-Spektroskopie).

WIP

