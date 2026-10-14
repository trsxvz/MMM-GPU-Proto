# Noyau perdu au premier lancement d'un noyau d'une bibliothèque partagée

ROCm 10.1 (HIP 7.16, paquets Python de TheRock), GPU intégré gfx1036
(Ryzen 7 9800X3D), MFEM 4.10 compilé avec HIP (scripts/build-hip-amd.sh).

Un noyau encore en attente sur le flux nul n'est jamais exécuté quand le
noyau suivant est le premier noyau de libmfem lancé par le programme, dont
le code est alors chargé (chargement différé, environ 2 ms). Le programme
annule des composantes d'un vecteur, comme `mfem::ParNonlinearForm::Mult`
pour les degrés de liberté imposés, puis lance `other = 3.` :

    make
    sg render -c ./lost-kernel            # 0 / 27 threads, y[0] = 1
    sg render -c "./lost-kernel warm"     # libmfem échauffée : 27 / 27
    sg render -c "AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3 ./lost-kernel"

Pas de perte avec un premier noyau de l'exécutable lui-même. Le chargement
immédiat (`HIP_ENABLE_DEFERRED_LOADING=0`) échoue : `libhipblaslt.so.1` n'a
pas de code pour gfx1036.

Effets dans mfem-mgis : résidu de Newton faux au premier pas (composantes
imposées non annulées, divergence), et fautes mémoire de hypre dans
rocsparse ou rocprim quand ce sont des copies ou des remplissages de hypre
qui sont perdus.

Avec la sérialisation, certains cas ne donnent plus un résultat faux mais
restent bloqués : dans les exemples rve et ssna303_3d_hypre sur le GPU, le
premier lancement d'un noyau de la bibliothèque des lois attend sans fin un
signal (`BusyWaitSignal::WaitAcquire` sous `hipLaunchKernel`, pile relevée
sous gdb).
