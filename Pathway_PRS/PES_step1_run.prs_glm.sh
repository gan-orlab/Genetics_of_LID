FILE=$1

Rscript ~/runs/yurils/lid-GWAS/PES/PRSice.R \
--prsice ~/runs/yurils/lid-GWAS/PES/PRSice_linux \
--base res.meta_no.${FILE}_glm.adj.tbl \
--binary-target T \
--gtf ~/runs/go_lab/gencode/gencode.v40lift37.annotation.gtf \
--thread 1 \
--target ${FILE} \
--beta \
--clump-p 0.05 \
--clump-kb 250 \
--pheno covar_glm.txt \
--pheno-col LID \
--cov covar_glm.txt \
--cov-col @PC[1-5],Sex,PD_AAO,LDOPA_cum,HY \
--prevalence 0.4 \
-p P-value \
--snp SNP \
--stat Effect \
--A1 Allele1 \
--A2 Allele2 \
--perm 10000 \
--msigdb /home/yurils/runs/yurils/lid-GWAS/PES/dopamine_pathways_and_LIDgenes.gmt \
--out ${FILE}_glm \
#--extract ${FILE}_glm.valid #uncomment this line if there are not valid SNP IDs
