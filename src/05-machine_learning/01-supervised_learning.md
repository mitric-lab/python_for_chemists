## Überwachtes Lernen

Liegen uns Datenpaare $\{(\vec{x}_i, y_i)\}_{i=1,2,\dots,N}$ vor, wobei $\vec{x}_i \in \mathbb{R}^n$ 
die Eingabe und $y_i \in \mathbb{R}$ die Ausgabe ist, so ist es nun unser Ziel, eine approximative
Abbildung $\hat{f}_{\theta} : \mathbb{R}^n \to \mathbb{R}$ zu finden, die die Eingabe auf die Ausgabe abbildet, 
d.h. $\hat{f}_{\theta}(\vec{x}_i) \approx y_i$ für alle $i$. Man bezeichnet $y_i$ auch als *Labels* oder *Targets*. Unser Modell $\hat{f}$ ist in 
der Regel parametrisiert, d.h. es gibt eine Menge von Parametern $\theta$, die die Abbildung definieren. 
Demnach ist unser Ziel, die Parameter $\theta$ so zu wählen, dass der Fehler zwischen der
tatsächlichen Ausgabe $y_i$ und der approximierten Ausgabe $\hat{f}_{\theta}(\vec{x}_i)$ minimiert wird.
Dieses Lernen wird als **überwachtes Lernen** bezeichnet, da wir $y_i$ für jede Eingabe 
$\vec{x}_i$ kennen und somit den Fehler direkt berechnen können.

Gemäß den obigen Überlegungen können wir festhalten, dass jedes ML-Modell, das auf überwachtem Lernen basiert,
die folgenden Komponenten enthält:

1. **Modell**: Die Funktion $\hat{f}_{\theta}$, die die Eingabe auf die Ausgabe abbildet.
2. **Verlustfunktion**: Die Funktion, die den Fehler zwischen der tatsächlichen Ausgabe und der approximierten Ausgabe misst.
3. **Optimierungsverfahren**: Der Algorithmus, der die Parameter $\theta$ des Modells so anpasst, dass der Fehler minimiert wird.

Dieses Setting lässt sich tatsächlich auch auf das unüberwachte Lernen übertragen, wobei die Verlustfunktion
durch eine allgemeine Objektivfunktion ersetzt wird, die es zu maximieren oder minimieren gilt.

### Regression

Für den Fall, dass die Labels $y_i$ kontinuierlich sind, also reelle Zahlen $y_i \in \mathbb{R}$ annnehmen, 
spricht man von **Regression**. Für den einfachsten Fall, die **lineare Regression**, und für $x_i \in \mathbb{R}$
haben