## Fourier-Transformation

Damit wir die Frequenzzerlegung für aperiodische Signale verstehen,
schauen wir uns die Fourier-Koeffizienten $c_n$ in 
Gl. {{eqref: eq:fourier_coeffs}} genauer an. Diese Koeffizienten geben an,
wie stark die Frequenz $n\omega$ in der Funktion $f(t)$ enthalten ist.
Wir können das als eine Abbildung verstehen:
$$
  c: \mathbb{Z} \ni n \mapsto c_n \in \mathbb{C}\,,
$$
oder explizit mit der Winkelfrequenz:
$$
  c: \mathbb{Z} \ni n \mapsto c(\omega_n) \in \mathbb{c}\,,
$$
wobei $\omega_n = n \omega$.

Eine aperiodische Funktion $f(t)$ kann als eine periodische Funktion
mit einer unendlich langen Periode betrachtet werden. Wir können also 
versuchen, $T \to \infty$ in Gl. {{eq: eq:fourier_coeffs}} zu setzen.
Wir erhalten
$$
  F_{\mathrm{b}}(\omega_n) = \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega_n t} \du t\,,
$$
wobei wir den Vorfaktor $1/T$ ignoriert haben, da dieser für $T \to \infty$
verschwinden würde und wir keine sinnvolle Definition für die Koeffizienten
haben würden. Deshalb werden diese Koeffizienten nicht mehr mit $c$, sondern
mit $F_{\mathrm{b}}$ bezeichnet. 

Die Grenzwertbildung hat auch eine Auswirkung auf die Winkelfrequenz 
$\omega_n$:
$$
  \omega_n = n \omega = \frac{2\pi n}{T} \to \omega\,.
$$
Durch $T \to \infty$ wird also der Abstand zwischen den Vielfachen
der Grundfrequenz infinitesimal klein, sodass wir von einer kontinuierlichen
Frequenzverteilung sprechen können. Deshalb verwenden wir die Notation
$\omega$ statt $\omega_n$. Das liefert uns die 
*Fourier-Transformation mit Rückwärtsnormierung*:
$$
  F_{\mathrm{b}}(\omega) = \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,.
  {{numeq}}{eq:fourier_transform_backward}
$$

Um die etwas komische Bezeichnung "Rückwärtsnormierung" zu verstehen, 
müssen wir zuerst die *Fourier-Synthese*, die Rekonstruktion der Funktion
aus den Frequenzanteilen, verstehen. Für eine periodische Ausgangsfunktion
ist die Fourier-Synthese durch die Fourier-Reihe in 
Gl. {{eq: eq:fourier_series}} gegeben. Da wir für eine aperiodische Funktion
ein kontinuierliches Frequenzspektrum haben, müssen wir die Summe durch ein
Integral ersetzen, also etwa wie
$$
  \mathring{f}(t) = \int_{-\infty}^{\infty} F_{\mathrm{b}}(\omega) \eu^{\iu \omega t} \du \omega\,,
$$
wobei $\mathring{f}(t)$ die rekonstruierte Funktion ist.
Setzen wir nun $F_{\mathrm{b}}(\omega)$ aus Gl. {{eq: eq:fourier_transform_backward}}
ein, erhalten wir
$$
  \begin{align}
    \mathring{f}(t) &= \int_{-\infty}^{\infty} \left( \int_{-\infty}^{\infty} f(t') \eu^{-\iu \omega t'} \du t' \right) \eu^{\iu \omega t} \du \omega \\
    &= \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
      f(t') \eu^{\iu \omega (t-t')} \du t' \du \omega \\
    &= \int_{-\infty}^{\infty} f(t') \left( \int_{-\infty}^{\infty} 
      \eu^{\iu \omega (t-t')} \du \omega \right) \du t'\,,
  \end{align}
$$
wobei wir unbekümmert die Reihenfolge der Integrationen vertauscht haben.

Das innere Integral sieht sehr ungewohnt aus, nach einer kurzen

```admonish derivation title="Rechnung für Interessierten" collapsible=true
Wir betrachten das Integral
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
$$
mit $\epsilon > 0$. Dieses Integral sieht erstmal noch komplizierter aus,
aber wir können die Exponentialfunktionen zusammenfassen und den Exponenten
quadratisch ergänzen:
$$
  \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} 
    = \eu^{\iu \omega t - \epsilon \omega^2} 
    = \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2 - t^2/4\epsilon}
    = \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2} \eu^{-t^2/4\epsilon}\,.
$$
Damit ist das obige Integral ein verschobenes Gauss-Integral:
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
  = \eu^{-t^2/4\epsilon} \int_{-\infty}^{\infty} \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2} \du \omega
  = \sqrt{\frac{\pi}{\epsilon}} \eu^{-t^2/4\epsilon} =: 2\pi d_{\epsilon}(t)\,.
$$

Da wir sowieso von $-\infty$ bis $\infty$ integrieren, spielt die Verschiebung
keine Rolle und wir erhalten das gleiche Ergebnis wie beim unverschobenen
Gauss-Integral. Das Ergebnis ist eine Funktion in Abhängigkeit von $t$,
die wir $2\pi d_{\epsilon}(t)$ nennen.

Dann schauen wir ein paar Eigenschaften der Funktion $d_{\epsilon}(t)$ an,
beginnend mit dem Integral:
$$
  \int_{-\infty}^{\infty} d_{\epsilon}(t) \du t
  = \frac{1}{2\pi}
  \int_{-\infty}^{\infty} \sqrt{\frac{\pi}{\epsilon}} \eu^{-t^2/4\epsilon} \du t
  = \frac{1}{2\pi} \sqrt{\frac{\pi}{\epsilon}} \sqrt{4\pi\epsilon}
  = 1\,.
$$
Das Integral dieser Funktion über die reellen Zahlen ist also $1$, unabhängig
von der Wahl vom $\epsilon$. 

Wie sieht es dann mit den Funktionswerten aus, wenn wir die Grenzwertbildung
$\epsilon \to 0$ durchführen? Die Funktion $d_{\epsilon}(t)$ ist eine
Gaußfunktion, welche in der $x$-Richtung mit $1/4\epsilon$ und in der
$y$-Richtung mit $\sqrt{\pi/\epsilon}$ skaliert ist. Für $\epsilon \to 0$
wird die Funktion also immer schmaler und höher, was zu der "Grenzfunktion"
$$
  d(t) := \lim_{\epsilon \to 0} d_{\epsilon}(t) = \begin{cases}
    \infty & x = 0 \\
    0 & x \neq 0
  \end{cases}
$$
führt. Obwohl der Funktionswert bei $0$ unendlich ist, haben wir vorher 
ausgerechnet, dass der Flächeninhalt unter der Kurve $1$ ist. Diese beiden
Eigenschaften definieren gerade die $\delta$-Funktion $delta{x}$, die unten 
im Ergebnis auftauchen. Genau genommen definieren diese zwei Eigenschaften 
keine Funktion, da es unklar ist, was $\delta(0)$ sein soll. Deshalb sollte 
man von einer $\delta$-Distribution sprechen. 

Die Grenzwertbildung $\epsilon \to 0$ führt auch dazu, dass der zusätzliche
Faktor im obigen Integral $\eu^{-\epsilon \omega^2}$ zu der konstanten 
Funktion $1$ wird. Damit erhalten wir das ursprüngliche Integral:
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \du \omega 
    = \lim_{\epsilon \to 0} \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
    = 2\pi \delta(t)\,.
$$
```

erhalten wir auch etwas "komisches":
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega (t-t')} \du \omega = 2\pi \delta(t-t')
$$
mit der 
[Dirac'schen Deltafunktion](https://de.wikipedia.org/wiki/Delta-Distribution)
$\delta(t)$.
Diese "Funktion" hat die Eigenschaft, dass das Integral über ihr Produkt
mit einer beliebigen Funktion $g(t)$ die ausgewertete Funktion an der Stelle
$0$ liefert:
$$
  \int_{-\infty}^{\infty} \delta(t) g(t) \du t = g(0)\,.
$$

Es gilt also insgesamt
$$
  \mathring{f}(t) = \int_{-\infty}^{\infty} f(t') (2\pi \delta(t-t')) \du t'
    = 2\pi f(t)\,.
$$
Die rekonstrurierte Funktion $\mathring{f}$ ist also um den Faktor $2\pi$
skaliert. Damit wir tatächlich die ursprüngliche Funktion rekonstruieren
können, müssen wir beim Rekonstruieren 
(rückwärts Fourier-Transformation oder inverse Fourier-Transformation)
die Funktion mit dem Normierungsfaktor $1/(2\pi)$ multiplizieren:
$$
  f(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F_{\mathrm{b}}(\omega) \eu^{\iu \omega t} \du \omega\,.
  {{numeq}}{eq:inverse_fourier_transform_backward}
$$
Diese Gleichung stellt nun die
*inverse Fourier-Transformation mit Rückwärtsnormierung* dar und
erklärt die Bezeichnung "Rückwärtsnormierung" in 
Gl. {{eq: eq:fourier_transform_backward}}, da die Normierung
bei der Rücktransformation berücksichtigt wird. Wir können gleichwertig
noch die *Fourier-Transformation mit Vorwärtsnormierung* definieren:
$$
  F_{\mathrm{f}}(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,,
  {{numeq}}{eq:fourier_transform_forward}
$$
wobei die Normierung bei der Vorwärts-Transformation berücksichtigt wird.
Es gibt auch noch die *Fourier-Transformation mit symmetrischer Normierung*:
$$
  F_{\mathrm{s}}(\omega) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,,
  {{numeq}}{eq:fourier_transform_symmetric}
$$
wobei der Normierungsfaktor bei Vorwärts- und Rückwärts-Transformation
aufgeteilt wird. 

Man findet in Lehrbüchern und wissenschaftlichen Arbeiten alle drei
Definitionen der Fourier-Transformation und diese sind alle gleichwertig.
Allerdings hat die symmetrische Normierung den Vorteil, dass die
Transformation eine
[unitäre Abbildung](https://de.wikipedia.org/wiki/Fourier-Transformation#Unitäre_Abbildung)
darstellt, welche Rechnungen in manchen Gebieten, wie z.B. in der 
Quantenmechanik, vereinfachen kann. 

Die Fourier-Transformation ist mathematisch schön und gut, aber ein Computer
kann mit einem kontinuierlichen Signal nicht viel anfangen, ganz zu schweigen
vom kontinuierlichen Frequenzspektrum. Deshalb führt uns zur
*diskrete Fourier-Transformation*, die für diskretisierte Signale
und Frequenzspektren geeignet ist. 

