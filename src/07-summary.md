# Zusammenfassung und Ausblick

In diesem Kurs haben Sie neben den Grundlagen des Programmierens mit Python 
auch einige grundlegende Konzepte und Methoden des maschinellen Lernens kennengelernt. 
Dazu zählen Methoden des überwachten und unüberwachten Lernens, sowie 
neuronale Netzwerke, wobei wir uns auf die Anwendung dieser Methoden in der Chemie 
konzentriert haben. Wir möchten dabei betonen, dass 
wir Ihnen in diesem Kurs nur einen kleinen Einblick in die Welt des 
maschinellen Lernens geben konnten. Mit Ihrem erworbenen Wissen und 
Fähigkeiten sind Sie jedoch nun durchaus in der Lage, auch komplexere 
Probleme zu lösen und eigene Projekte zu realisieren, wozu wir Sie 
ausdrücklich ermutigen. 

In der Chemie spielen maschinelle Lernverfahren eine immer 
wichtigere Rolle, da sie es ermöglichen, komplexe Zusammenhänge in großen 
Datenmengen zu erkennen und zu nutzen. Wir möchten Ihnen daher abschließend 
einen Einblick in den aktuellen Stand des maschinellen Lernens in der Chemie 
geben, sowie einige aktuelle Forschungsthemen aus unserem Arbeitskreis vorstellen.

## Aktueller Stand des maschinellen Lernens in der Chemie

Um eines der spannendsten Anwendungsgebiete des maschinellen Lernens (in der Chemie 
und allgemein), das *generative Modellieren*, zu motivieren, möchten wir 
Sie zunächst . 
Wir haben gesehen, dass neuronale Netzwerke sowohl zur Klassifikation als auch 
zur Regression eingesetzt werden können. Dabei benötigen wir zum Training 
des Netzwerks 

In der Chemie finden sie daher 
Anwendung in der Vorhersage von Moleküleigenschaften, wie z.B. der Synthetisierbarkeit 
oder der Aktivität von Wirkstoffen. 

Können wir neuronale Netzwerke jedoch auch 
für das unüberwachte Lernen einsetzen? 

## Aktuelle Forschungsthemen in unserem Arbeitskreis

Im Folgenden möchten wir Ihnen einen Einblick in die aktuellen Forschungsthemen 
am Lehrstuhl für Theoretische Chemie geben. Dazu zählen die Arbeitskreise von 
Prof. Dr. Roland Mitrić und Dr. Merle Röhr. Zudem stellen wir Ihnen kurze 
Codebeispiele vor (Python und andere Programmiersprachen), die beispielhaft 
für die Forschungsthemen stehen und in denen 
Sie sicherlich einige der in diesem Kurs erlernten Konzepte wiedererkennen werden.

### Nicht-adiabadische Dynamik mit Trajectory Surface Hopping

Im Rahmen der Trajectory Surface Hopping Methode bewegen sich die Moleküle 
in einem $3N$-dimensionalen Konfigurationsraum, wobei $N$ die Anzahl der Atome
im Molekül ist. Selbst bei kleinen Molekülen ist dieser Konfigurationsraum
für Menschen nur schwer vorstellbar. Allerdings ist ein Großteil dieses 
Konfigurationsraums für die Dynamik
uninteressant, da sie wegen ihrer hohen Energie (z.B. durch sehr große
Bindungsabsände) für das Molekül nicht zugänglich sind. Das erlaubt
uns, Dimensionalitätsreduktionstechniken zu verwenden, um den 
hochdimensionalen Konfigurationsraum auf eine niedrigdimensionale
Darstellung zu projizieren, ohne dabei zu viele Informationen zu verlieren.
Im folgenden Codebeispiel wird die multidimensionale Skalierung (MDS)
verwendet, um eine niedrigdimensionale Darstellung des Konfigurationsraums
zu berechnen. 

```python
def perform_traj_mds(atnos, ref_coords, traj_coords, ndim=2, 
                     npi_pairs=None, dist_pairs=None, max_dist=None):
    nref = len(ref_coords)
    coords = np.concatenate((ref_coords, traj_coords))
    descriptors, weights = get_geom_descriptor(
        atnos, coords, 
        npi_pairs=NPI_PAIRS, dist_pairs=DIST_PAIRS, max_dist=MAX_DIST,
    )
    dissimilarities = get_desc_distmat(descriptors, weights)

    mds = MDS(
        n_components=ndim, dissimilarity='precomputed', random_state=42,
        n_init=32, max_iter=1000, eps=1e-4,
    )
    embedding = mds.fit_transform(dissimilarities)

    return embedding
```

Diese Technik wurde für die Untersuchung der Photodissoziationsdynamik
von Cyclobutanon eingesetzt.[^miao2024] Eine Darstellung von 
4 repräsentativen Trajektorien sind in der folgenden Abbildung zu sehen.

![Cyclobutanon](./assets/figures/07-summary/selected_mds.svg)

### Fragment-basierte Methoden für die Berechnung von angeregten Zuständen in großen molekularen Aggregaten

Die theoretische Beschreibung von dynamischen Prozessen, wie Exciton-Transfer 
oder Ladungstransfer, in organischen Halbleitern und komplexen molekularen 
System erfordert Methoden, die die Berechnung von großen molekularen 
Aggregaten ermöglichen, die aus tausenden Atomen bestehen. Die 
"traditionellen" quantenchemische Methoden, wie z.B. Hartree-Fock oder
Dichtefunktionaltheorie, sind aufgrund ihrer starken Skalierung mit der
Systemgröße nicht in der Lage, solche Systeme in angemessener Zeit zu
berechnen. Um diese Problematik anzugehen, entwickelt unsere Gruppe 
neue theoretische Methoden, die solche großen molekularen Systeme beschreiben 
können. Dazu kombinieren wir semiempirische quantenchemische Methoden (DFTB) 
mit einem Fragmentierungsansatz (FMO) in einem neuen theoretischen 
Formalismus, der es erlaubt die angeregten Zustände von großen molekularen 
Aggregaten zu berechnen und die Molekulardynamik dieser angeregten Zustände 
zu simulieren.[^einsele2023]<sup>,</sup>[^einsele2024]

```rust
// create the A matrix from the orbital energy differences, 
// the Coulomb and the exchange contributions
let h: Array2<f64> = self.fock_and_coulomb() - self.exchange();
// solve the eigenvalue problem A x = w A using the eigenvalue decomposition
let (eigenvalues, eigenvectors) = h.eigh(UPLO::Upper).unwrap();

// Reference to the o-v transition charges.
let q_ov: ArrayView2<f64> = self.properties.q_ov().unwrap();

// The transition charges for all excited states are computed.
let q_trans: Array2<f64> = q_ov.dot(&eigenvectors);

// The Mulliken transition dipole moments are computed.
let tr_dipoles: Array2<f64> = mulliken_dipoles(q_trans.view(), &self.atoms);

// The oscillator strengths are computed.
let f: Array1<f64> = oscillator_strength(eigenvalues.view(), tr_dipoles.view());
```

![FMO](./assets/figures/07-summary/fig_10.svg)

### Simulation von stark gekoppelten Licht-Materie-Systemen

Energietransport in excitonischen Materialien spielt eine große Rolle für die 
Anwendung in vielen optoelektronischen Systemen. Während in der Vergangenheit 
versucht wurden, den Energietransport durch strukturelle Veränderung der verwendeten 
Moleküle zu verbessern, zielt ein Forschungsschwerpunkt unseres Arbeitskreises 
darauf ab, dies durch Kopplung der 
elektronischen Übergänge an starke elektromagnetische Felder zu erreichen, bspw. 
in Mikrokavitäten. Die Quasiteilchen, die in solchen Systemen entstehen, werden 
Polaritonen genannt. Genauer beschäftigen wir uns mit der theoretischen Beschreibung
von Polaritonen. Dazu müssen viele Konzepte aus dem Studium der theoretischen 
Chemie zum Einsatz gebracht werden und zusätzlich kombiniert werden mit Methoden 
der Quantenelektrodynamik. Um reale Systeme zu berechnen, werden die entstehenden 
Gleichungen numerisch gelöst.

```python
system.build_system()
system.set_operators()
system.set_H() 
e, v = np.linalg.eigh(np.real(system.H))
coeff = np.dot(v.T, np.dot(system.a + system.a_dagger, v))
```

### Theoretische Untersuchung von kleinen Metallclustern

Dieser Forschungsbereich konzentriert sich auf kleine Metallcluster 
sowohl im Zeit- als 
auch im Energieraum (zeitabhängige und zeitunabhängige Prozesse). Ziel ist es, diese 
kleinen Metallcluster theoretisch so genau wie möglich zu beschreiben. Obwohl diese 
Cluster nur eine geringe Anzahl von Atomen enthalten, ist ihre Untersuchung aufgrund 
der komplizierten elektronischen Natur der d- und f-Schalen der Metalle sehr komplex. 
Diese Cluster sind besonders interessant, weil das Verständnis ihrer katalytischen 
Aktivität stark von theoretischen Studien profitieren könnte, die derzeit nicht in 
ausreichendem Maße verfügbar sind.

Dieses Code-Snippet stammt aus einem Programm, das die eindimensionale Schrödinger-Gleichung 
numerisch exakt löst. Diese Gleichung ist fundamental in der Quantenmechanik, insbesondere 
für das Verständnis von zweiatomigen Molekülen. 
 
```python
def getHamiltonian(self):
    # Initialisiert eine Hamilton-Matrix mit komplexen Nullen
    self.hamiltonian = np.zeros((dim, dim), dtype=complex)
    
    # Berechnet die Impulswerte pk
    pk = (2. * np.pi / (dim * self.deltax)) * (indrange - dim / 2.)
    
    # Berechnet den kinetischen Term tk
    tk = (pk**2) / (2.0 * self.m)
    
    # Berechnet den Exponentialterm W
    W = np.exp(2 * np.pi * 1.0j * indrange / dim)
    
    # Schleife durch jede Zeile i der Hamilton-Matrix
    for i in range(dim):
        # Berechnet die temporären Werte für die Fourier-Transformation
        tmp = tk * (W**i) * ((-1)**i)
        
        # Führt die Fourier-Transformation durch und aktualisiert die Hamilton-Matrix-Zeile
        self.hamiltonian[i,:] = oneminusone * np.fft.fft(tmp) / dim
```

Der vorliegende Code hat die Aufgabe, eine Hamilton-Matrix (`self.hamiltonian`)
zu konstruieren und eine Fourier-Transformation des kinetischen Teils 
durchzuführen, da dieser Operator im Impulsraum multiplikativ ist. Zunächst 
wird die Dimension des Grids (`dim`) ermittelt und eine komplexe Nullmatrix 
für den Hamiltonian (`self.hamiltonian`) initialisiert. Anschließend werden 
die Impulswerte (`pk`) berechnet und daraus der kinetische Term (`tk`) abgeleitet. 
Durch eine exponentielle Funktion (`W`) und eine Schleife über die Matrixzeilen 
wird der kinetische Term in den Impulsraum transformiert und mittels 
Fourier-Transformation (`np.fft.fft`) in die Hamilton-Matrix eingetragen. 


![Cer](./assets/figures/07-summary/ce2_together.svg)
Abbildung: Theoretische Untersuchung von Ce<sub>2</sub>. 
Links: Vergleich zwischen dem beobachteten NeNePo-Signal im Experiment 
und den theoretischen Simulationen (Zeitraum).
Mitte: Berechnete potentielle Energiekurven. Die hervorgehobenen Zustände 
stellen die identifizierten Zustände dar, die für das beobachtete Signal 
verantwortlich sind.
Rechts: Vergleich zwischen dem experimentellen Photoelektronenspektrum von 
Ce<sub>2</sub> und dem simulierten Spektrum (Energieraum).


### Optimierung von Dimeren für Singlet Fission

Singlet Fission ist ein Prozess in bestimmten organischen Molekülen, bei dem ein einzelnes 
angeregtes Singulett-Exciton in zwei Triplett-Excitonen zerfällt. Dies kann theoretisch den 
Wirkungsgrad von Solarzellen erhöhen, da ein einzelnes Photon zwei Elektronen-Loch-Paare 
anregen kann.

```python
scaler = StandardScaler()
X = scaler.fit_transform(X)
pca = PCA(n_components=n_components)
pca.fit(X)
transformed_data = pca.transform(X)
explained_variance_ratios = np.array(pca.explained_variance_ratio_)
# Get the principal components
components = np.array(pca.components_)
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(transformed_data)
labels = np.array(kmeans.labels_)
centroids = np.array(kmeans.cluster_centers_)
```

Wir haben zufällige Dimere konstruiert und diese hinsichtlicher einer Eigenschaft, der 
sogenannten Singlet-Fission Rate optimiert und die optimierten Dimer-Strukturen, dann 
analysiert wie z.B Translation und Rotation usw.. Um dann in den Daten Gemeinsamkeiten 
zu finden, habe ich eine PCA und ein K-Means-Clustering genutzt.

![PCA](./assets/figures/07-summary/3D_cluster.svg)

### Supramolekulare Aggregate des Lichtsammelkomplexes der grünen Schwefelbakterien

Mein Forschungsbereich konzentriert sich auf die Strukturaufklärung von Aggregaten, die in 
grünen Schwefelbakterien vorkommen. Diese Bakterien sind bemerkenswert für ihre Fähigkeit, 
in extremen Lichtverhältnissen Photosynthese zu betreiben, wobei sie spezielle 
Lichtsammelkomplexe, sogenannte Chlorosomen, nutzen. Diese Chlorosomen bestehen aus dicht 
gepackten Bacteriochlorophyll-Molekülen, die außergewöhnlich effiziente Energietransferprozesse 
aufweisen. Diese Aggregate bestehen oft aus mehreren konzentrischen Ringen, die eine 
zylindrische Form zeigen (siehe Abbildung).

![Chlorosomen](./assets/figures/07-summary/SimplifiedPaperTube_combined_wholeSide_gray_cropped.svg)

Der Code Ausschnitt zeigt eine Implementierung des Frenkel-Exziton-Hamiltonian, der die 
Wechselwirkungen zwischen den Übergangsdipolmomenten der einzelenen Moleküle im Aggregat 
einbezieht. Damit kann ein Spektrum des Aggregates simuliert werden und dieses mit experimentellen 
Daten verglichen werden.

```python
def getExcitonHamiltonian(self):
        nmol = len(self.allMolecules)
        self.hamiltonian = np.zeros((nmol, nmol))
        for i, m in enumerate(self.allMolecules):
            self.hamiltonian[i, i] = m.siteEnergy
        # convert units! a.u/Angstrom --> eV/a.u
        self.hamiltonian += self.getExcitonCoupling() * toang * toev

    # create exciton states for aggregate
    def getExcitonStates(self):
        print("First step: Calculate hamiltonian.")
        self.getExcitonHamiltonian()
        print("Second step: Solve hamiltonian.")
        self.excitonStateEnergies, self.excitonStateCoefficients = np.linalg.eigh(self.hamiltonian)
        self.getSiteDipoles()
        self.tdMoments = []
        print("Third step: tdMoments.")
        self.tdMoments = np.dot(np.transpose(self.excitonStateCoefficients), self.dipoles)
        np.save("td", self.tdMoments)
        np.save("excitonStates", self.excitonStateEnergies)
```


---



[^miao2024]: X. Miao, K. Diemer, R. Mitrić, *J. Chem. Phys.* **2024**, *160*, 124309.

[^einsele2023]: R. Einsele, J. Hoche, R. Mitrić, *J. Chem. Phys.* **2023**, *158*, 044121.

[^einsele2024]: R. Einsele, R. Mitrić, *J. Comput. Theor. Chem.* **2024**, just accepted.
