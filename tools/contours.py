import cv2
import numpy as np

# 1. Charger l'image PNG (sans canal alpha)
img = cv2.imread("sim_HOLO_PHASE.png")  # Format BGR standard

# Vérifier si l'image est chargée
if img is None:
    print("Erreur : impossible de charger l'image")
    exit()

# 2. Créer une image en niveaux de gris pour la détection de contours
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# 3. Améliorer le contraste (optionnel mais recommandé)
# Appliquer une égalisation d'histogramme
gray = cv2.equalizeHist(gray)

# Ou appliquer un seuillage adaptatif si l'image est peu contrastée
# _, gray = cv2.threshold(gray, 127, 255, cv2.THRESH_BINARY)

# 4. Détection des bords avec Canny
# Ajuster les seuils selon votre image (essayez 50/150 ou 30/100 si nécessaire)
edges = cv2.Canny(gray, 20, 100)

# 5. Optionnel : dilater les contours pour les rendre plus visibles
kernel = np.ones((3, 3), np.uint8)
edges = cv2.dilate(edges, kernel, iterations=1)

# 6. Trouver les contours
contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

print(f"Nombre de contours détectés : {len(contours)}")

# 7. Dessiner les contours sur une copie de l'image originale
output = img.copy()

# Dessiner en couleur vive (rouge épais) pour bien voir
cv2.drawContours(output, contours, -1, (0, 0, 255), thickness=3)  # Rouge, épaisseur 3

# 8. Sauvegarder les résultats
cv2.imwrite("sim_HOLO_PHASE_contours.png", output)
cv2.imwrite("sim_HOLO_PHASE_edges.png", edges)  # Sauvegarder aussi les bords détectés

print("Images sauvegardées avec succès")

# 9. Optionnel : afficher les images pour vérifier
# cv2.imshow("Original", img)
# cv2.imshow("Edges", edges)
# cv2.imshow("Contours", output)
# cv2.waitKey(0)
# cv2.destroyAllWindows()

import matplotlib.pyplot as plt

plt.figure(figsize=(15, 5))

plt.subplot(2, 2, 1)
plt.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
plt.title("Original")
plt.axis('off')

plt.subplot(2, 2, 2)
plt.imshow(gray, cmap='gray')
plt.title("Niveaux de gris")
plt.axis('off')

plt.subplot(2, 2, 3)
plt.imshow(edges, cmap='gray')
plt.title("Edges")
plt.axis('off')

plt.subplot(2, 2, 4)
plt.imshow(cv2.cvtColor(output, cv2.COLOR_BGR2RGB))
plt.title("Contours")
plt.axis('off')

plt.tight_layout()
plt.show()
