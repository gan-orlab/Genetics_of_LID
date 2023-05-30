FILE=$1

Rscript ~/runs/yurils/lid-GWAS/PES/PRSice.R \
--prsice ~/runs/yurils/lid-GWAS/PES/PRSice_linux \
--base res.meta_no.${FILE}_cox.adj.tbl \
--gtf ~/runs/go_lab/gencode/gencode.v40lift37.annotation.gtf \
--thread 1 \
--target ${FILE} \
--beta \
--clump-p 0.05 \
--clump-kb 250 \
--pheno covar_cox.txt \
--pheno-col LID \
--cov covar_cox.txt \
--cov-col @PC[1-5],Sex,PD_AAO,LDOPA_total,HY \
-p P-value \
--snp SNP \
--stat Effect \
--A1 Allele1 \
--A2 Allele2 \
--perm 10000 \
--msigdb /home/yurils/runs/yurils/lid-GWAS/PES/dopamine_pathways_and_LIDgenes.gmt \
--out ${FILE}_cox \
#--extract ${FILE}_cox.valid  #uncomment this line if there are not valid SNP IDs
