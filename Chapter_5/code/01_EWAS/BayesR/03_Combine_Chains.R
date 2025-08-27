# Take last 250 rows from each chain, keeping header from one file
  traits=(bmi body_fat Glucose HDL_cholesterol Total_cholesterol whr)
  for trait in ${traits[@]}
  do   
head -n+1 Sigma/${trait}_chain1.csv > Sigma/head_sigma.csv
head -n+1 Beta/${trait}_chain1.csv > Beta/head_beta.csv
head -n+1 Comp/${trait}_chain1.csv > Comp/head_comp.csv
for chain in {1..4}
do
tail -250 Sigma/${trait}_chain${chain}.csv >> Sigma/${trait}_temp.csv
tail -250 Beta/${trait}_chain${chain}.csv >> Beta/${trait}_temp.csv
tail -250 Comp/${trait}_chain${chain}.csv >> Comp/${trait}_temp.csv
done
cat Sigma/head_sigma.csv Sigma/${trait}_temp.csv > Sigma/${trait}_processed.csv
cat Beta/head_beta.csv Beta/${trait}_temp.csv > Beta/${trait}_processed.csv
cat Comp/head_comp.csv Comp/${trait}_temp.csv > Comp/${trait}_processed.csv
rm */*temp.csv
rm */head*.csv
done

### R ####
setwd("")
pdf("convergence_fullmodel.pdf",width=12,height=11)
par(mfrow=c(4,6))
for(trait in c("bmi","body_fat","Glucose","HDL_cholesterol","Total_cholesterol","whr")){
  for(chain in 1:4){
    tmp = read.csv(paste0(trait, "_chain", chain, ".csv"))
    plot(rowSums(tmp), ylab=paste0("Sigma Chain ", chain), main = paste0(trait, " Chain ", chain))
  }
}
dev.off()
