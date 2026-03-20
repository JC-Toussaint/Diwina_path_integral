#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colormaps
from PIL import Image
from PIL.ExifTags import TAGS

from matplotlib.colors import ListedColormap
from PyQt5.QtWidgets import (
    QApplication, QGraphicsView, QGraphicsScene, QGraphicsPixmapItem,
    QMainWindow, QWidget, QVBoxLayout, QLabel, QHBoxLayout, QTabWidget,
    QGridLayout
)
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen, QFont
from PyQt5.QtCore import Qt, QRectF, pyqtSignal, QObject

class ZoomSynchronizer(QObject):
    """Gestionnaire centralisé pour synchroniser les zooms entre toutes les vues"""
    zoom_changed = pyqtSignal(float, object)  # factor, sender
    
    def __init__(self):
        super().__init__()
        self.viewers = []
        self._updating = False  # Flag pour éviter les boucles infinies
    
    def register_viewer(self, viewer):
        """Enregistre un viewer pour la synchronisation"""
        self.viewers.append(viewer)
        viewer.zoom_applied.connect(self.on_zoom_changed)
    
    def on_zoom_changed(self, factor, sender):
        """Propage le zoom à tous les autres viewers"""
        if self._updating:
            return
            
        self._updating = True
        try:
            for viewer in self.viewers:
                if viewer != sender:
                    viewer.apply_zoom_from_sync(factor)
        finally:
            self._updating = False

class AxisWidget(QLabel):
    def __init__(self, orientation='horizontal', parent=None):
        super().__init__(parent)
        self.orientation = orientation
        self.min_coord = 0
        self.max_coord = 100
        if orientation == 'horizontal':
            self.setFixedHeight(20)
        else:
            self.setFixedWidth(45)
        
    def update_range(self, min_coord, max_coord):
        self.min_coord = min_coord
        self.max_coord = max_coord
        self.update()
        
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Configuration du pinceau et de la police
        pen = QPen(Qt.black, 1)
        painter.setPen(pen)
        font = QFont("Arial", 12)
        painter.setFont(font)
        
        rect = self.rect()
        
        if self.orientation == 'horizontal':
            # Axe horizontal
            y_pos = rect.height() - 1
            
            # Ligne principale
            painter.drawLine(0, y_pos, rect.width(), y_pos)
            
            # Calculer les ticks
            range_coord = self.max_coord - self.min_coord
            if range_coord > 0:
                # Nombre fixe de ticks : 11
                n_ticks = 11
                tick_values = np.linspace(self.min_coord, self.max_coord, n_ticks)
                
                for i, value in enumerate(tick_values):
                    x_pos = int(i * rect.width() / (n_ticks - 1))
                    
                    # Dessiner le tick vers le haut
                    painter.drawLine(x_pos, y_pos, x_pos, y_pos - 3)
                    
                    # Dessiner le label au-dessus
                    text = f"{int(value)}"
                    text_rect = painter.fontMetrics().boundingRect(text)
                    painter.drawText(x_pos - text_rect.width() // 2, 
                                   y_pos - 5, text)
        else:
            # Axe vertical
            x_pos = rect.width() - 1
            
            # Ligne principale
            painter.drawLine(x_pos, 0, x_pos, rect.height())
            
            # Calculer les ticks
            range_coord = self.max_coord - self.min_coord
            if range_coord > 0:
                # Nombre fixe de ticks : 11
                n_ticks = 11
                tick_values = np.linspace(self.min_coord, self.max_coord, n_ticks)
                
                for i, value in enumerate(tick_values):
                    y_pos = int(i * rect.height() / (n_ticks - 1))
                    
                    # Dessiner le tick vers la gauche
                    painter.drawLine(x_pos, y_pos, x_pos - 3, y_pos)
                    
                    # Dessiner le label à gauche
                    text = f"{int(value)}"
                    text_rect = painter.fontMetrics().boundingRect(text)
                    painter.drawText(x_pos - text_rect.width() - 5, 
                                   y_pos + text_rect.height() // 2, text)

class ImageViewer(QGraphicsView):
    # Signal émis quand la vue change
    view_changed = pyqtSignal(float, float, float, float)  # min_x, max_x, min_y, max_y
    # Signal émis quand un zoom est appliqué
    zoom_applied = pyqtSignal(float, object)  # factor, sender
    
    def __init__(self, image_path):
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.pixmap = QPixmap(image_path)
        self.pixmap_item = QGraphicsPixmapItem(self.pixmap)
        self.scene().addItem(self.pixmap_item)

        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.setSceneRect(QRectF(self.pixmap.rect()))
        self.fitInView(self.sceneRect(), Qt.KeepAspectRatio)
        
        # Variables pour la synchronisation
        self._sync_updating = False
        
        # Émettre le signal initial
        self.emit_view_changed()

#    def showEvent(self, event):
#        super().showEvent(event)
#        self.fitInView(self.sceneRect(), Qt.KeepAspectRatio)
#        self.emit_view_changed()
    
    def wheelEvent(self, event):
        if self._sync_updating:
            return
            
        zoom_in_factor = 1.01
        zoom_out_factor = 1 / zoom_in_factor
        factor = zoom_in_factor if event.angleDelta().y() > 0 else zoom_out_factor
        
        self.scale(factor, factor)
        self.emit_view_changed()
        
        # Émettre le signal de zoom pour synchronisation
        self.zoom_applied.emit(factor, self)

    def apply_zoom_from_sync(self, factor):
        """Applique un zoom depuis la synchronisation (sans émettre de signal)"""
        self._sync_updating = True
        try:
            self.scale(factor, factor)
            self.emit_view_changed()
        finally:
            self._sync_updating = False

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pos = self.mapToScene(event.pos())
#            print(f"Clicked at scene position: {pos.x():.1f}, {pos.y():.1f}")
        super().mousePressEvent(event)

    def scrollContentsBy(self, dx, dy):
        super().scrollContentsBy(dx, dy)
        self.emit_view_changed()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.emit_view_changed()

    def emit_view_changed(self):
        # Obtenir le rectangle visible dans le viewport
        viewport_rect = self.viewport().rect()
        
        # Mapper les coins du viewport vers les coordonnées de la scène
        top_left = self.mapToScene(viewport_rect.topLeft())
        top_right = self.mapToScene(viewport_rect.topRight())
        bottom_left = self.mapToScene(viewport_rect.bottomLeft())
        bottom_right = self.mapToScene(viewport_rect.bottomRight())
        
        # Obtenir les limites de la zone visible
        min_x = min(top_left.x(), bottom_left.x())
        max_x = max(top_right.x(), bottom_right.x())
        min_y = min(top_left.y(), top_right.y())
        max_y = max(bottom_left.y(), bottom_right.y())
        
        # Limiter aux dimensions de l'image
        img_rect = self.pixmap.rect()
        min_x = max(0, min_x)
        max_x = min(img_rect.width(), max_x)
        min_y = max(0, min_y)
        max_y = min(img_rect.height(), max_y)
        
        self.view_changed.emit(min_x, max_x, min_y, max_y)

def read_image_metadata(image_path):
    """
    Lit les métadonnées d'une image PNG et extrait les valeurs min/max si disponibles
    """
    try:
        with Image.open(image_path) as img:
            metadata = img.text if hasattr(img, 'text') else {}
            
            # Chercher les valeurs dans les métadonnées textuelles
            vmin, vmax, unit = None, None, None
            
            for key, value in metadata.items():
                if 'min' in key.lower():
                    try:
                        vmin = float(value)
                    except:
                        # Essayer d'extraire un nombre de la chaîne
                        import re
                        match = re.search(r'[-+]?\d*\.?\d+', value)
                        if match:
                            vmin = float(match.group())
                elif 'max' in key.lower():
                    try:
                        vmax = float(value)
                    except:
                        import re
                        match = re.search(r'[-+]?\d*\.?\d+', value)
                        if match:
                            vmax = float(match.group())
                elif 'unit' in key.lower():
                    unit = value
                        
#            print(f"Métadonnées pour {image_path}:")
#            for key, value in metadata.items():
#                print(f"  {key}: {value}")
                
            return vmin, vmax, unit
            
    except Exception as e:
        print(f"Error while reading metadata for {image_path}: {e}")
        return None, None, None

def read_custom_colormap(path):
    data = np.loadtxt(path)
    x = data[:, 0]  # Positions x
    rgb = data[:, 1:4]  # Composantes RGB
    
    # Normaliser les positions x entre 0 et 1
    x_norm = (x - x.min()) / (x.max() - x.min())
    
    # Créer des points d'interpolation pour une colormap lisse
    n_colors = 256
    x_interp = np.linspace(0, 1, n_colors)
    
    # Interpoler chaque composante RGB
    r_interp = np.interp(x_interp, x_norm, rgb[:, 0])
    g_interp = np.interp(x_interp, x_norm, rgb[:, 1])
    b_interp = np.interp(x_interp, x_norm, rgb[:, 2])
    
    # Combiner les composantes RGB interpolées
    rgb_interp = np.column_stack([r_interp, g_interp, b_interp])
    
    # S'assurer que les valeurs sont dans [0, 1]
    rgb_interp = np.clip(rgb_interp, 0, 1)
    
    return ListedColormap(rgb_interp)

def generate_colorbar_pixmap(cmap, vmin=0, vmax=1, unit=None):
    # Créer une figure avec des dimensions appropriées
    fig, ax = plt.subplots(figsize=(1.2, 4), dpi=100)
    
    # Supprimer le contenu de l'axe principal
    ax.remove()
    
    # Créer un axe pour la colorbar avec plus d'espace pour les labels
    cbar_ax = fig.add_axes([0.1, 0.1, 0.2, 0.8])  # [left, bottom, width, height]
    
    # Créer la colorbar avec les vraies valeurs
    # La colormap sera interpolée linéairement entre vmin et vmax
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
    cb.outline.set_visible(False)
    
    # Formater les labels avec 2 chiffres après la virgule
    cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2e}'))
    
    # Ajuster la taille de police et la position des labels
    cb.ax.tick_params(labelsize=8)
    cb.ax.yaxis.set_label_position("right")
    cb.ax.yaxis.tick_right()
    
    # Définir des ticks explicites pour une meilleure répartition
    n_ticks = 6
    tick_values = np.linspace(vmin, vmax, n_ticks)
    cb.set_ticks(tick_values)
    
    # Ajouter l'unité si disponible
    if unit:
        cb.set_label(unit, rotation=90, labelpad=0)
    
    # Ajuster les marges pour laisser de la place aux labels
    fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05)
    
    # Convertir la figure en image
    fig.canvas.draw()
    
    # Obtenir les données de l'image
    buf = fig.canvas.buffer_rgba()
    w, h = fig.canvas.get_width_height()
    img = np.frombuffer(buf, dtype=np.uint8).reshape((h, w, 4))
    
    # Convertir RGBA en RGB
    img_rgb = img[:, :, :3].copy()
    
    plt.close(fig)

    # Créer QImage et QPixmap - convertir en bytes
    img_bytes = img_rgb.tobytes()
    image = QImage(img_bytes, w, h, w * 3, QImage.Format_RGB888)
    return QPixmap.fromImage(image)

class ImageTab(QWidget):
    def __init__(self, image_path, cmap=None, vmin=None, vmax=None, unit=None):
        super().__init__()
        
        # Layout principal en grille
        layout = QGridLayout()
        
        # Créer les axes
        self.top_axis = AxisWidget('horizontal')
        self.left_axis = AxisWidget('vertical')
        
        # Créer le visualiseur d'image
        self.viewer = ImageViewer(image_path)
        
        # Connecter le signal de changement de vue aux axes
        self.viewer.view_changed.connect(self.update_axes)
        
        # Créer la colorbar si une colormap est fournie
        colorbar_widget = None
        if cmap:
            # Utiliser les valeurs par défaut si pas de métadonnées
            if vmin is None:
                vmin = 0
            if vmax is None:
                vmax = 1
                
            colorbar_pixmap = generate_colorbar_pixmap(cmap, vmin, vmax, unit)
            colorbar_widget = QLabel()
            colorbar_widget.setPixmap(colorbar_pixmap)
            colorbar_widget.setScaledContents(False)
            colorbar_widget.setAlignment(Qt.AlignCenter)
            colorbar_widget.setMaximumWidth(100)
        
        # Organiser les widgets dans la grille
        # Position: (row, col, rowspan, colspan)
        layout.addWidget(self.top_axis, 0, 1, 1, 1)  # Axe horizontal en haut
        layout.addWidget(self.left_axis, 1, 0, 1, 1)  # Axe vertical à gauche
        layout.addWidget(self.viewer, 1, 1, 1, 1)  # Visualiseur d'image au centre
        if colorbar_widget:
            layout.addWidget(colorbar_widget, 1, 2, 1, 1)  # Colorbar à droite
        
        # Ajuster les proportions des colonnes et lignes
        layout.setColumnStretch(0, 0)  # Colonne de l'axe vertical : taille fixe
        layout.setColumnStretch(1, 1)  # La colonne de l'image s'étend
        if colorbar_widget:
            layout.setColumnStretch(2, 0)  # Colonne colorbar : taille fixe
        
        layout.setRowStretch(0, 0)     # Ligne de l'axe horizontal : taille fixe  
        layout.setRowStretch(1, 1)     # La ligne de l'image s'étend
        
        # Réduire les marges et l'espacement pour maximiser l'espace image
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        
        self.setLayout(layout)
    
    def update_axes(self, min_x, max_x, min_y, max_y):
        """Met à jour les axes avec les nouvelles coordonnées visibles"""
        self.top_axis.update_range(min_x, max_x)
        self.left_axis.update_range(min_y, max_y)

class MainWindow(QMainWindow):
    def __init__(self, image_paths, colorbar_cmap_list, metadata_list):
        super().__init__()
        self.setWindowTitle("pathIntegral png viewer")

        # Créer le synchroniseur de zoom
        self.zoom_sync = ZoomSynchronizer()

        # Créer le widget central avec des onglets
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        # Créer un onglet pour chaque image
        for idx, (path, cmap) in enumerate(zip(image_paths, colorbar_cmap_list)):
            # Récupérer les métadonnées pour cette image
            vmin, vmax, unit = metadata_list[idx] if idx < len(metadata_list) else (None, None, None)
            
            # Créer l'onglet avec l'image et sa colorbar
            tab = ImageTab(path, cmap, vmin, vmax, unit)
            
            # Enregistrer le viewer pour la synchronisation du zoom
            self.zoom_sync.register_viewer(tab.viewer)
            
            # Extraire le nom du fichier pour le titre de l'onglet
            tab_name = path.split("/")[-1].replace(".png", "").replace("sim_", "")
            
            # Ajouter l'onglet au widget d'onglets
            self.tab_widget.addTab(tab, tab_name)

        # Redimensionner la fenêtre
        # Redimensionner la fenêtre et empêcher le redimensionnement
        self.resize(1000, 900)
        self.setFixedSize(1000, 900)
        
if __name__ == "__main__":
    app = QApplication(sys.argv)

    image_list = [
        "sim_HOLO_PHASE.png",
        "sim_MZ.png",
        "sim_PATH_LENGTH.png",
        "sim_STXM_XMCD.png"
    ]

    gray_map = colormaps.get_cmap("gray")
    colorbar_cmap_list = [
        gray_map,
        gray_map,
        gray_map,
        gray_map
    ]

    # Lire les métadonnées de toutes les images
    metadata_list = []
    for image_path in image_list:
        vmin, vmax, unit = read_image_metadata(image_path)
        metadata_list.append((vmin, vmax, unit))
#        print(f"{image_path} - Min: {vmin}, Max: {vmax}, Unit: {unit}")

    window = MainWindow(image_list, colorbar_cmap_list, metadata_list)
    window.show()
    sys.exit(app.exec_())
