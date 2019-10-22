#a partir du dossier Seconde_partie ./../Troisieme_partie/07_test_cluster_GO.sh > ../Troisieme_partie/resultats_04_comp_meth.txt
i=1
while read -r line
do
echo cluster $i hypergeometrique
./03_blastset.py --query "$line" --sets /home/guest/Desktop/Cheryn_ALI_Projet_IHD/Premiere_partie/RESULTATS/GOterm_prot_dir_Solanum_tuberosum.sets --adjust --alpha 0.05 --measure 'hypergeometric'
echo cluster $i binomial
./03_blastset.py --query "$line" --sets /home/guest/Desktop/Cheryn_ALI_Projet_IHD/Premiere_partie/RESULTATS/GOterm_prot_dir_Solanum_tuberosum.sets --adjust --alpha 0.05 --measure 'binomial'
./03_blastset.py --query "$line" --sets /home/guest/Desktop/Cheryn_ALI_Projet_IHD/Premiere_partie/RESULTATS/GOterm_prot_dir_Solanum_tuberosum.sets --adjust --alpha 0.05 --measure 'chi2'
i=$(($i+1))
done < /home/guest/Desktop/Cours_M2/IDH/Projet/troixieme_partie/RESULTATS/GN_clusters.txt