#/usr/bin/bash
workplace=$(pwd)
mkdir -p ${workplace}/Gymnosperm_odb10
##After filtering and  spliting low-copy orthologs
ls | grep "OG[0-9]*.fa" |sed 's/.fa//' > OG_list
#Gene names in BUSCO v5.4.7 should start with number, rename OGs
for n in $(cat ${workplace}/OG_list);do
	newname=$(echo ${n} | sed 's/OG[0]*//')
	mv ${n}.fa ${newname}.fa
	cat ${newname}.fa >> ${workplace}/Gymnosperm_odb10/refseq_db.faa
done
## generating the required profiles
mkdir -p ${workplace}/Gymnosperm_odb10/hmms
mkdir -p ${workplace}/Gymnosperm_odb10/prfl
for i in $(cat ${workplace}/OG_list | sed 's/OG[0]*//');do
	mafft --auto ${i}.fa > ${i}.mafft
	hmmbuild ${i}.hmm ${i}.mafft
    hmmsearch ${i}.hmm ${i}.mafft > ${i}.hmm.searchout 
	cat ${i}.hmm.searchout | sed -n '15,25p' > ${i}.hmm.search.stats
	hmmemit ${i}.hmm > ${i}.ancestral6  
	hmmemit -N 10 ${i}.hmm > ${i}.ancestral_variants
	${Augustus_scripts}/msa2prfl.pl ${i}.mafft --qij=${Augustus_profile}/default.qij > ${workplace}/Gymnosperm_odb10/prfl/${i}.prfl
	mv ${i}.hmm ${workplace}/Gymnosperm_odb10/hmms
	sed 's/-.*//' ${i}.ancestral >> ${workplace}/Gymnosperm_odb10/ancestral
	sed 's/-sample/_/g' ${i}.ancestral_variants >> ${workplace}/Gymnosperm_odb10/ancestral_variants
done
## Score and length cut off
mkdir -p ${workplace}/score
mkdir -p ${workplace}/length
for m in $(cat ${workplace}/OG_list | sed 's/OG[0]*//');do
	cat ${m}.hmm.search.stats | awk '{print $2}' | grep '^[0-9]' > ${workplace}/score/${m}.score
	##the 90% value of the lowest score and two significant digits
	a=$(cat ${workplace}/score/${m}.score | sort -k1n | head -n 1); b=$(echo "$a*0.9"|bc); echo $(awk 'BEGIN{printf "%.2f",'$b'}') > ${workplace}/score/${m}.0.9.lowest.score
	sed "s/^/${m}\t/g" ${workplace}/score/${m}.0.9.lowest.score >> ${workplace}/Gymnosperm_odb10/scores_cutoff
	##the length range was set to the mean and twice the standard deviation
	c=$(seqkit fx2tab ${m}.fa -l -n |awk '{print $2}' |awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}');d=$(echo "$c*2"|bc); echo $(awk 'BEGIN{printf "%.2f",'$d'}') >> ${workplace}/length/stdevp_length
	e=$(seqkit fx2tab ${m}.fa -l -n |awk '{print $2}' | awk '{total+=$1} END {print total/NR}') ; echo $(awk 'BEGIN{printf "%.2f",'$e'}') >>${workplace}/length/mean_length
	echo -e "${m}\t0" >>${workplace}/length/part_length
done
paste -d "\t" ${workplace}/length/part_length ${workplace}/length/stdevp_length ${workplace}/length/mean_length > ${workplace}/Gymnosperm_odb10/lengths_cutoff
## set dataset.cfg
sed -i 's/embryophyta_odb10$/gymnosperm_odb10/' dataset.cfg
sed -i 's/arabidopsis$/gymnosperm/' dataset.cfg
sed -i 's/1614$/1603/' dataset.cfg
sed -i 's/50$/7/' dataset.cfg
sed -i 's/90000$/1600000/' dataset.cfg
sed -i 's/120000$/2200000/' dataset.cfg
## Sampling for assessing the reliability
## Genome and protein data of Gnetum montanum can be found at https://doi.org/10.5061/dryad.ht76hdrdr
## Cite The Welwitschia genome reveals a unique biology underpinning extreme longevity in deserts
## Randomly removing BUSCOs from Gymnosperm(gym) and embryophyta(emb) geneset
## seed for seqkit sample
function rand(){
 min=$1
 max=$(($2-$min+1))
 num=$(date +%s%N)
 echo $(($num%$max+$min))
}
for i in $(seq 0.1 0.2 0.5);do
	mkdir -p ${workplace}/BUSCO${i}
	cd ${workplace}/BUSCO${i}
	for n in $(seq 1 1 5);do
			rnd=$(rand 1 100)
			seqkit sample -n ${i} -s ${rnd} ${workplace}/Gn.gym.fa |grep ">" |sed 's/>//' >Gn.gym_${i}_${n}_list
			seqkit sample -n ${i} -s ${rnd} ${workplace}/Gn.emb.fa |grep ">" |sed 's/>//' >Gn.emb_${i}_${n}_list
			done
			##remove samping BUSCOS sequence in pep
			seqkit grep -f Gn.gym_${i}_${n}_list -v ${workplace}/Gnetum.pep.fa > Gnetum.pep.gym.${i}_${n}.fa
			seqkit grep -f Gn.emb_${i}_${n}_list -v ${workplace}/Gnetum.pep.fa > Gnetum.pep.emb.${i}_${n}.fa
			##mask samping BUSCOS sequence in genome 
			for m in $(cat Gn.gym_${i}_${n}_list);do
				grep "${m};" ${workplace}/G.montanum.final.gff >>Gn.gym_${i}_${n}.gff
			done
			for o in $(cat Gn.emb_${i}_${n}_list);do
				grep "${o};" ${workplace}/G.montanum.final.gff >>Gn.emb_${i}_${n}.gff
			done
			bedtools maskfasta -fi ${workplace}/G.montanum.final.fasta -bed Gn.gym_${i}_${n}.gff -fo G.montanum.final.gym.${i}.${n}.fasta
			bedtools maskfasta -fi ${workplace}/G.montanum.final.fasta -bed Gn.emb_${i}_${n}.gff -fo G.montanum.final.emb.${i}.${n}.fasta	
done
##Randomly removing sequences from protein data of Gnetum montanum
for i in $(seq 0.1 0.2 0.5);do
	mkdir -p ${workplace}/Protein${i}
	cd ${workplace}/Protein${i}
		for n in $(seq 1 1 5);do
			rnd=$(rand 1 100)
			seqkit sample -n ${i} -s ${rnd} ${workplace}/Gnetum.pep.fa | grep ">" | sed 's/>//' >Gn.pep_${i}_${n}_list
			##remove samping pep sequence
			seqkit grep -f Gn.pep_${i}_${n}_list -v ${workplace}/Gnetum.pep.fa > Gnetum.pep.${i}_${n}.fa
			##mask samping pep sequence in genome 
				for m in $(cat Gn.pep_${i}_${n}_list);do
					grep "${m};" ${workplace}/G.montanum.final.gff >>Gn.pep_${i}_${n}.gff
				done
			bedtools maskfasta -fi ${workplace}/G.montanum.final.fasta -bed Gn.pep_${i}_${n}.gff -fo G.montanum.final.pep.${i}.${n}.fasta
		done
done
