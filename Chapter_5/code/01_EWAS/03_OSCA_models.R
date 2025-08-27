########### OSCA models 20k####################

## example OSCA script  

## loop through phenotypes 
for i in *.phen
do  
# extracts phenotype name 
A=$( echo $i | cut -d"." -f1)


#runs OSCA
osca_Linux \
--linear \
--befile osca_20k \
--pheno $i \
--qcovar covariate.cov \
--fast-linear \
--out /results/${A}_full \
--methylation-m
done 

