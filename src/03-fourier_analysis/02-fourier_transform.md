## Fourier-Transformation

Bevor wir auf die Frequenzzerlegung für aperiodische Signale eingehen,
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
mit einer unendlich langen Periode betrachtet werden. Unser Ansatz ist demnach, 
in Gl. {{eqref: eq:fourier_coeffs}} $T \to \infty$ zu setzen.
Wir erhalten
$$
  F_{\mathrm{b}}(\omega_n) = \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega_n t} \du t\,,
$$
wobei wir den Vorfaktor $1/T$ ignoriert haben, da dieser für $T \to \infty$
verschwindet. Da wir so keine sinnvolle Definition für die Koeffizienten erhalten,
werden diese nicht mehr als $c$, sondern als $F_{\mathrm{b}}$ bezeichnet. 

Die Grenzwertbildung hat auch eine Auswirkung auf die Winkelfrequenz 
$\omega_n$:
$$
  \omega_n = n \omega = \frac{2\pi n}{T} \to \omega\,.
$$
Durch $T \to \infty$ wird der Abstand zwischen den Vielfachen
der Grundfrequenz infinitesimal klein, sodass wir von einer kontinuierlichen
Frequenzverteilung sprechen können. Deshalb verwenden wir die Notation
$\omega$ statt $\omega_n$. Das liefert uns die 
*Fourier-Transformation mit Rückwärtsnormierung*:
$$
  F_{\mathrm{b}}(\omega) = \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,.
  {{numeq}}{eq:fourier_transform_backward}
$$

Um die wenig intuitive Bezeichnung "Rückwärtsnormierung" zu verstehen, 
müssen wir zuerst die *Fourier-Synthese*, die Rekonstruktion der Funktion
aus den Frequenzanteilen, diskutieren. Für eine periodische Ausgangsfunktion
ist die Fourier-Synthese durch die Fourier-Reihe in 
Gl. {{eqref: eq:fourier_series}} gegeben. Da wir für eine aperiodische Funktion
ein kontinuierliches Frequenzspektrum haben, müssen wir die Summe durch ein
Integral ersetzen, also etwa wie
$$
  \mathring{f}(t) = \int_{-\infty}^{\infty} F_{\mathrm{b}}(\omega) \eu^{\iu \omega t} \du \omega\,,
$$
wobei $\mathring{f}(t)$ die rekonstruierte Funktion ist.
Setzen wir nun $F_{\mathrm{b}}(\omega)$ aus Gl. {{eqref: eq:fourier_transform_backward}}
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

Das innere Integral sieht sehr ungewohnt aus. Nach einer kurzen Rechnung
erhalten wir
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega (t-t')} \du \omega = 2\pi \delta(t-t')
$$
mit der [Dirac'schen Deltafunktion](https://de.wikipedia.org/wiki/Delta-Distribution) 
$\delta(t)$.

```admonish derivation title="Rechnung für Interessierten" collapsible=true
Wir betrachten das Integral
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
$$
mit $\epsilon > 0$. Wir können die Exponentialfunktionen zusammenfassen und den Exponenten
quadratisch ergänzen:
$$
  \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} 
    = \eu^{\iu \omega t - \epsilon \omega^2} 
    = \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2 - t^2/4\epsilon}
    = \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2} \eu^{-t^2/4\epsilon}\,.
$$
Damit ist das obige Integral ein verschobenes Gauss-Integral, dessen analytische 
Lösung bekannt ist:
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
  = \eu^{-t^2/4\epsilon} \int_{-\infty}^{\infty} \eu^{-\epsilon (\omega - \iu t/2\epsilon)^2} \du \omega
  = \sqrt{\frac{\pi}{\epsilon}} \eu^{-t^2/4\epsilon} =: 2\pi d_{\epsilon}(t)\,.
$$

Da wir sowieso von $-\infty$ bis $\infty$ integrieren, spielt die Verschiebung
keine Rolle und wir erhalten das gleiche Ergebnis wie beim unverschobenen
Gauss-Integral. Das Ergebnis ist eine Funktion in Abhängigkeit von $t$,
die wir $2\pi d_{\epsilon}(t)$ nennen.

Im Folgenden betrachten wir einige Eigenschaften der Funktion $d_{\epsilon}(t)$, 
beginnend mit dem Integral:
$$
  \int_{-\infty}^{\infty} d_{\epsilon}(t) \du t
  = \frac{1}{2\pi}
  \int_{-\infty}^{\infty} \sqrt{\frac{\pi}{\epsilon}} \eu^{-t^2/4\epsilon} \du t
  = \frac{1}{2\pi} \sqrt{\frac{\pi}{\epsilon}} \sqrt{4\pi\epsilon}
  = 1\,.
$$
Das Integral dieser Funktion über die reellen Zahlen ist also normiert, unabhängig
von der Wahl von $\epsilon$. 

Wie sieht es dann mit den Funktionswerten $d_{\epsilon}(t)$ aus, wenn wir die 
Grenzwertbildung $\epsilon \to 0$ durchführen? Die Funktion $d_{\epsilon}(t)$ ist eine
Gaußfunktion, dessen Definitionsbereich mit $1/4\epsilon$ und dessen Zielmenge 
mit $\sqrt{\pi/\epsilon}$ skaliert ist. Für $\epsilon \to 0$
wird die Funktion also immer schmaler und höher, was zu der "Grenzfunktion"
$$
  d(t) := \lim_{\epsilon \to 0} d_{\epsilon}(t) = \begin{cases}
    \infty & x = 0 \\
    0 & x \neq 0
  \end{cases}
$$
führt. Obwohl der Funktionswert bei $0$ unendlich ist, haben wir vorher 
berechnet, dass der Flächeninhalt unter der Kurve $1$ ist. Diese beiden
Eigenschaften definieren gerade die $\delta$-Funktion $\delta{x}$, die wir im
Folgenden verwenden werden. Genau genommen definieren diese beiden Eigenschaften 
keine Funktion, da unklar bleibt, was $\delta(0)$ bedeutet. Deshalb sollte 
man eher von einer $\delta$-Distribution sprechen. 

Die Grenzwertbildung $\epsilon \to 0$ führt ebenfalls dazu, dass der zusätzliche
Faktor $\eu^{-\epsilon \omega^2}$ im obigen Integral zu $1$ wird. Damit erhalten 
wir:
$$
  \int_{-\infty}^{\infty} \eu^{\iu \omega t} \du \omega 
    = \lim_{\epsilon \to 0} \int_{-\infty}^{\infty} \eu^{\iu \omega t} \eu^{-\epsilon \omega^2} \du \omega
    = 2\pi \delta(t)\,.
$$
```

Diese "Funktion" hat die Eigenschaft, dass das Integral über ihr Produkt
mit einer beliebigen Funktion $g(t)$, diese Funktion ausgewertet an der Stelle
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
können, müssen wir während dieser Rekonstruktion 
(rückwärts Fourier-Transformation oder inverse Fourier-Transformation)
die Funktion mit dem Normierungsfaktor $1/(2\pi)$ multiplizieren:
$$
  f(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F_{\mathrm{b}}(\omega) \eu^{\iu \omega t} \du \omega\,.
  {{numeq}}{eq:inverse_fourier_transform_backward}
$$
Diese Gleichung stellt nun die
*inverse Fourier-Transformation mit Rückwärtsnormierung* dar und
erklärt die Bezeichnung der 
Gl. {{eqref: eq:fourier_transform_backward}}, da die Normierung
bei der Rücktransformation berücksichtigt wird. Wir können analog
die *Fourier-Transformation mit Vorwärtsnormierung* definieren:
$$
  F_{\mathrm{f}}(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,,
  {{numeq}}{eq:fourier_transform_forward}
$$
wobei die Normierung bei der Vorwärts-Transformation berücksichtigt wird. 
Die entsprechende inverse Transformation enthält demnach
keinen Normierungsfaktor.
Es gibt auch noch die *Fourier-Transformation mit symmetrischer Normierung*:
$$
  F_{\mathrm{s}}(\omega) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} f(t) \eu^{-\iu \omega t} \du t\,,
  {{numeq}}{eq:fourier_transform_symmetric}
$$
wobei der Normierungsfaktor auf die Vorwärts- und Rückwärts-Transformation
aufgeteilt wird. 

Man findet in Lehrbüchern und wissenschaftlichen Arbeiten alle drei
Definitionen der Fourier-Transformation. Alle diese Definitionen sind gleichwertig
und korrekt, sofern die zugehörige inverse Transformation entsprechend definiert wird.
Allerdings hat die symmetrische Normierung den Vorteil, dass die
Transformation eine
[unitäre Abbildung](https://de.wikipedia.org/wiki/Fourier-Transformation#Unitäre_Abbildung)
darstellt, was Rechnungen in manchen Gebieten, wie z.B. in der 
Quantenmechanik, vereinfachen kann. 

Die Fourier-Transformation ist aus mathematischer Sicht schön und gut, aber ein Computer
kann mit einem kontinuierlichen Signal nicht viel anfangen, ganz zu schweigen
von kontinuierlichen Frequenzspektren. Deshalb führt uns der nächste Abschnitt zur
*diskreten Fourier-Transformation*, die wir in der Praxis verwenden können.

