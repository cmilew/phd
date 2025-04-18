import numpy as np


def calculate_dta_gamma_index(coords1, coords2, dta_threshold):
    """
    Compare deux séries de coordonnées avec un critère de Distance To Agreement (DTA).

    Args:
        coords1 (array-like): Liste de tuples (x, y) pour la première série de coordonnées.
        coords2 (array-like): Liste de tuples (x, y) pour la deuxième série de coordonnées.
        dta_threshold (float): Le seuil de distance pour le critère DTA.

    Returns:
        float: Le gamma index moyen pour les deux séries de coordonnées.
    """
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # Initialiser une liste pour stocker les valeurs du gamma index
    gamma_indexes = []

    # Comparer chaque point de coords1 avec chaque point de coords2
    for point1 in coords1:
        min_distance = np.inf
        for point2 in coords2:
            # Calculer la distance euclidienne entre point1 et point2
            distance = np.linalg.norm(point1 - point2)
            # Trouver la distance minimale
            if distance < min_distance:
                min_distance = distance
        # Calculer le gamma index pour ce point
        gamma_index = min_distance / dta_threshold
        gamma_indexes.append(gamma_index)

    # Calculer la valeur moyenne du gamma index
    mean_gamma_index = np.mean(gamma_indexes)

    return mean_gamma_index


# Exemple d'utilisation
coords1 = [(1, 1), (2, 2), (3, 3)]
coords2 = [(1.1, 1.1), (2.1, 2.1), (3.1, 3.1)]
dta_threshold = 0.5

gamma_index = calculate_dta_gamma_index(coords1, coords2, dta_threshold)
print(f"Gamma Index moyen: {gamma_index}")
