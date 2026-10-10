# Newton inexact sur la poutre, jacobiens par noyau (10/10/2026)

Mêmes cas que `results/L5-poutre-2026-10-10-0429/` (poutre r4, 111 843
inconnues, 10 pas, tous les noyaux, script `cas/mesures-L5-ew.sh`, tolérance
initiale d'Eisenstat-Walker 0,5), avec le calcul des jacobiens par un noyau
(mfem-mgis 8b0b47acfb ; `resume.txt` indique le commit de base). Charge 1,1
au départ.

| Cas | GPU, avant | GPU | CPU 8 processus, avant | CPU 8 processus |
|---|---|---|---|---|
| Newton | 3,12 s | 2,98 s | 13,59 s | 13,57 s |
| Newton, Eisenstat-Walker 1 | 3,01 s | 2,74 s | 15,73 s | 15,76 s |
| Newton, Eisenstat-Walker 2 | 3,32 s | 3,10 s | 17,70 s | 17,82 s |
| NewtonK0 | 2,57 s | 2,56 s | 11,56 s | 11,59 s |
| NewtonK0, Eisenstat-Walker 1 | 2,14 s | **1,90 s** | 10,93 s | **11,00 s** |
| NewtonK0, Eisenstat-Walker 2 | 2,19 s | 1,99 s | 12,30 s | 12,26 s |

Sur le GPU, le gain va de 0 à 0,27 s selon le cas, une exécution par cas :
la géométrie de 197 000 tétraèdres devait coûter environ 0,2 s sur l'hôte
(un huitième des 1,8 s mesurées en r5). Sur 8 processus, chaque processus
n'en calcule qu'un huitième : pas de différence mesurable.
