# Titre du graphique
set title "Temps d'exécution en fonction de la taille de la matrice"

# Libellés des axes
set xlabel "Taille de la matrice"
set ylabel "Temps d'exécution (s)"

# Tracé de la courbe à partir du fichier
plot "exc.txt" using 1:2 with linespoints title 'Algorithme ALPHA'
