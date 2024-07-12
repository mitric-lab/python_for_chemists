# Zusammenfassung und Ausblick

## Aktueller Stand des maschinellen Lernens in der Chemie

## Aktuelle Forschungsthemen in unserem Arbeitskreis


### Nicht-adiabadische Dynamik mit Trajectory Surface Hopping

Im Rahmen der Trajectory Surface Hopping Methode bewegen sich die Moleküle 
in einem $3N$-dimensionalen Konfigurationsraum, wobei $N$ die Anzahl der Atome
im Molekül ist. Selbst bei kleinen Molekülen ist dieser Konfigurationsraum
für Menschen nur schwer vorstellbar.

Allerdings ist ein Großteil dieses Konfigurationsraums für die Dynamik
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

---











[^miao2024]: X. Miao, K. Diemer, R. Mitrić, *J. Chem. Phys.* **2024**, *160*, 124309.
