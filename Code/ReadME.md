# Projet LookBackOption 
  - Étudier et implémenter la méthode de calcul du delta (et du gamma) par calcul de
Malliavin proposée par [BGK03] pour différentes options de type lookback de votre
choix.
   ▷ Calculer également le delta par la méthode du flot, par différences finies : à pas fixe, à
   pas décroissant. Dans le second cas on choisira avec soin la constante "c" qui gouverne
   le pas de différence finie. On pourra ajouter des réducteurs de variance.
   ▷ Comparer le poids de Malliavin avec celui obtenu par l’approche martingale dans [Gob04].
   ▷ Calculer le gamma par différences finies, par une approche mixte (différences finies sur
   la méthode du flot) et comparer les résultats obtenus par les deux méthodes. On propose
   d’ajouter une méthode de réduction de de variance

## Nécessaire pour la compilation
 - Il faut installer cmake
 - Il faut installer la librairie armadillo
 - Il faut installer python, numpy et matplotlib

## Compilation
 - cmake CMakeLists.txt
 - make

## Exécution
 - Exécuter en lançant le fichier ./LookBackOption 252 20000
 - 252 = nombre de pas 
 - 20000 = nombre de simulations
