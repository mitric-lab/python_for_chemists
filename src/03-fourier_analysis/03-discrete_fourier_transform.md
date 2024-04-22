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
Dadurch erhalten wir die *diskrete Fourier-Transformation* (DFT).

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
Ohne genau auf die Theorie der FTIR-Spektroskopie einzugehen, gilt
$$
  S(\widetilde{\nu}) \propto \mathcal{F}\{I(x)\}(\widetilde{\nu})\,,
$$
wobei $I(x)$ das gemessene Interferogramm und $S(\tilde{\nu})$ das gewünschte
IR-Spektrum ist. Das Symbol $\mathcal{F}$ steht für die 
(diskrete oder kontinuierliche) Fourier-Transformation.
In anderen Worten, wenn wir das gemessene Interferogramm Fourier-transformieren,
erhalten wir das IR-Spektrum, da uns der Vorfaktor nicht interessiert.

Ihr findet hier das Interferogramm der
<a href="../codes/03-fourier_analysis/ir_bg.txt" download>Hintergrundmessung</a> und der
<a href="../codes/03-fourier_analysis/ir_spl.txt" download>Probemessung</a>.
Der Anfang der Datei der Hintergrundmessung (`ir_bg.txt`) sieht wie folgt aus:
```
{{#include ../codes/03-fourier_analysis/ir_bg.txt:0:10}}
```
Die Datei enthält zwei Spalten, wobei die erste die Verschiebung in Anzahl der
Schritten und die zweite die Intensität ist. Die Datei der Probemessung
(`ir_spl.txt`) sieht ähnlich aus. Bei dem verwendeten Spektrometer beträgt
die Schrittweite 0.3165&nbsp;&mu;m.

Als erstes importieren wir die benötigten Bibliotheken und definieren die
Schrittweite als Konstante `XSTEP`:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:imports}}
```
Danach importieren wir die Messdaten:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:import_data}}
```
Da wir später die Hintergrundmessung von der Probemessung abziehen wollen,
sollen wir sicherstellen, dass die Verschiebungen in beiden Messungen
übereinstimmen, weshalb wir es mit dem `assert`-Befehl überprüft haben.
Wir können die Inteferogramme anschließend plotten:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:plot_interferograms}}
```
Um die Details zu erkennen, haben wir den Plotbereich auf die Schritte von
$-250$ bis $250$ beschränkt. Das Diagramm sehen wir hier:
![Interferogramme der Hintergrund- und Probemessung](../assets/figures/03-fourier_analysis/ir_interferograms.svg)
Beide Interferogramme sehen extrem ähnlich aus.

Um das IR-Spektrum zu erhalten, müssen wir das Differenzsignal 
Fourier-transformieren:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:dft_signal}}
```
Nachdem wir das Differenzsignal berechnet haben, wurde das Nullarray
`int_nu` mit der Länge der Messdaten erstellt. Beachte dabei, dass wir den
Datentyp `complex` verwenden müssen, da die Fourier-Transformation
komplexe Zahlen zurückgibt. Danach haben wir mit 
[np.arange](https://numpy.org/doc/stable/reference/generated/numpy.arange.html)
das Array $(0, 1, \cdots, N-1)$ erstellt, welches wir genutzt haben, um nach
Gl. {{eqref: eq:dft_backward}} durch eine Schleife über `k` die
diskrete Fourier-Transformation zu berechnen.

Danach berechnen wir die "Frequenzen" nach Gl. {{eqref: eq:dft_frequencies}}.
Hier sollten wir aber eine kleine Anpassung vornehmen: Fa unser Ausgangssignal
in Verschiebung $x$ gemessen wurde und nicht in der Zeit $t$, ist die 
reziproke Größe die Wellenzahl $\widetilde{\nu} = 1/x$ bzw. die Länge des
Wellenvektors $k = 2\pi/x$. Weil das IR-Spektrum konventionell in Abhängigkeit
der Wellenzahl dargestellt wird, sollten wir den Vorfaktor $2\pi$ in
Gl. {{eqref: eq:dft_frequencies}} durch $1$ ersetzen. Das führt zu der 
folgenden Implementierung:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:dft_freq}}
```
Hier haben wir mit der `if`-`else`-Anweisung den Fallunterschied zwischen
geradem und ungeradem `nx` berücksichtigt. Dabei wurden den positiven
und negativen Frequenzanteilen separat berechnet. Am Ende haben wir die
mit `np.concatenate` zusammengefügt.

Weil es angenehmer ist, mit monotonen Wellenzahlen zu arbeiten, sortieren
wir im Anschluss der Fourier-Transformation die Wellenzahlen und das Spektrum:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:dft_shift}}
```
Die Funktion
[np.argsort](https://numpy.org/doc/stable/reference/generated/numpy.argsort.html)
gibt die Indizes der sortierten Werte zurück. Mit diesen Indizes können wir
dann alle Arrays gemäß der Sortierung eines Arrays (hier `nu`) sortieren.

Wenn wir uns an "Herleitung" der DFT erinnern, haben wir dabei angenommen, 
dass die Funktion $f(t)$ an den Stellen $t = 0, 1, \cdots, N-1$ gesampelt wurde.
In unserem Fall startet die Verschiebung aber bei $-8191 \cdot 0.3165$&nbsp;\[&mu;m\]
und das Intervall zwischen den Schritten beträgt 0.3165&nbsp;\[&mu;m\].
Deshalb müssen wir die ursprüngliche Funktion skalieren und verschieben, um
die korrekte Fourier-Transformation zu erhalten. Dabei helfen uns die
folgenden Beziehungen:
$$
  \begin{align}
    \mathcal{F}\{f(t - b)\} &= \eu^{-\iu b \omega} F(\omega)\,,
    {{numeq}}{eq:fourier_properties_shift} \\
    \mathcal{F}\{f(a t)\} &= \frac{1}{|a|} F\left(\frac{\omega}{a}\right)\,,
    {{numeq}}{eq:fourier_properties_scale}
  \end{align}
$$
wobei $F(\omega) = \mathcal{F}\{f(t)\}$. Diese Beziehungen werden Sie in der
Übung zeigen. Die Gl. {{eqref: eq:fourier_properties_scale}} 
besagt, dass eine Skalierung des Arguments der Funktion $f(t)$ zu einer
Skalierung des Arguments und des Betrags der Fourier-Transformierten führt.
Während diese Beziehung den Vorfaktor $1/\Delta t$ in Gl. {{eqref: eq:dft_frequencies}}
erklärt, können wir sie in diesem Fall für die DFT vom Signal ignorieren, da uns
die absolute Skalierung des Spektrums nicht interessiert. 
Die Gl. {{eqref: eq:fourier_properties_shift}} besagt, dass eine Verschiebung
im Zeitraum um $b$ zu einer Multiplikation der Fourier-Transformierten mit
dem Phasenfaktor $\eu^{-\iu b \omega}$ führt. Weil die DFT die 
Ausgangsfunktion ab $0$ auswertet, müssten wir sie um 
$-8191 \cdot 0.3165$&nbsp;\[&mu;m\] verschieben, was in unserer 
Implementierung in `x_grid[0]` gespeichert ist. Das führt zu der folgenden
Manipulation im Fourier-Raum
```python
{{#include ../codes/03-fourier_analysis/ftir.py:signal_shift}}
```
Hier haben wir den *In-place*-Operator `*=` verwendet, um die Werte von
`int_nu` direkt zu ändern. Es kann gezeigt werden, dass die 
Fourier-Transformierte einer symmetrischen Funktion reell ist, was wir 
hier mit dem `assert`-Befehl überprüfen. Dabei wird das Attribut `imag`
eines komplexen Arrays verwendet, um den Imaginärteil zu erhalten.

Eine letzte kosmetische Anpassung ist die Änderung der Einheit der Variable
`nu` von $\mathrm{m^{-1}}$ zu $\mathrm{cm^{-1}}$:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:nu_conversion}}
```
Hier wurde der In-place-Operator `/=` verwendet, um die Werte von `nu` direkt
durch $100$ zu teilen.

Endlich können wir das IR-Spektrum plotten:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:plot_spectrum}}
```
Wir haben mit dem Attribute `real` den Realteil des Spektrums erhalten.
Außerdem haben wir den Plotbereich auf $0$ bis $4000\ \mathrm{cm^{-1}}$ beschränkt.
Das Spektrum sieht wie folgt aus:
![IR-Spektrum der Probe](../assets/figures/03-fourier_analysis/ir_spectrum.svg)

Wir erkennen deutlich die aromatische C-C-Streckschwingung bei etwa 
$1600\ \mathrm{cm^{-1}}$ sowie aromatische und aliphatische 
C-H-Streckschwingungen um $3000\ \mathrm{cm^{-1}}$. Weiterhin sehen wir
noch Gerüstschwingungen im Fingerprint-Bereich. Es handelt sich hierbei
um ein IR-Spektrum von Mesitylen.

Bei einem so verbreiteten Algorithmus wie die diskrete Fourier-Transformation
bietet NumPy selbstverständlich eine Implementierung. Diese können wir
mit der Funktion
[`np.fft.fft`](https://numpy.org/doc/stable/reference/generated/numpy.fft.fft.html)
wie folgt aufrufen:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:numpy_fft}}
```
Das Umsortieren der Frequenzen und des Spektrums erfolgt mit der Funktion
[`np.fft.fftshift`](https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html):
```python
{{#include ../codes/03-fourier_analysis/ftir.py:numpy_fftshift}}
```

Es wird Ihnen vielleicht auffallen, dass die NumPy-Implementierung viel 
schneller ist als unsere eigene Implementierung. Das liegt einerseits daran,
dass NumPy in C/C++ geschrieben ist und deshalb schneller ist als Python. 
Andererseits liegt es daran, dass in NumPy eine effiziente Implementierung der
diskreten Fourier-Transformation, die sog. *Fast Fourier Transformation* (FFT), 
verwendet wird. Unsere Implementierung der DFT hat eine Komplexität von
$\mathcal{O}(N^2)$, d.h., dass für große Anzahlen von Datenpunkten $N$ die
die Laufzeit quadratisch mit $N$ steigt. Die FFT hat dagegen eine Komplexität
von $\mathcal{O}(N \log N)$.

Die FFT ist keine Nährung der DFT, sondern eine exakte Implementierung. Wir 
können die gleichen Skalierungen an den Intensitäten und Frequenzen der FFT
anwenden und die Ergebnisse der unserer Implementierung vergleichen:
```python
{{#include ../codes/03-fourier_analysis/ftir.py:numpy_verification}}
```

