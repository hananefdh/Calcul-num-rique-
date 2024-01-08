# Ouverture du fichier de données
set datafile separator "\n"

# Configuration du tracé
set terminal png  # Format de sortie du graphique
set output "plot.png"  # Nom du fichier de sortie

# Configuration du titre et des libellés des axes
set title "historique convergence RESVEC"
set xlabel "itérations"
set ylabel "data"

# Tracer les données
plot 'plot_data.txt' using 1 with lines title 'Data'
