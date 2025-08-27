## Load required libraries (may have to install first)

library(qqman)
library(data.table)


# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis

#read in results file
bmi = fread("bmi_EWAS_full_model_filtered_CpGs.csv") 
bmi_ewas = as.data.frame(bmi)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 
## saving out a tiff image for the qq plot
tiff("BMI_Full_model_QQ_Plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## qq plot 
qq(bmi_ewas$p, main = expression("Q-Q plot of BMI (kg/m"^2*") OSCA EWAS"),  pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis

text(norm_x,norm_y, "\u03BB = 6.7")
## save out - should say 'null device' once done, if not, run dev.off() once more  and check file is there 
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 
# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-bmi_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)

lambda
#[1] 6.74566

## Step 3 - make Manhattan plot 
ppi = 300
tiff("/BMI_Full_model_Manhat_Plot.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  bmi_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D", "#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = expression("BMI(kg/m"^2*")")
)

dev.off()

#########################################################################################################################################################################
#########################################################################################################################################################################
## Read in EWAS File 

whr = fread("whr_EWAS_full_model_filtered_CpGs.csv") 
whr_ewas = as.data.frame(whr)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 

## start saving out a tiff image for the qq plot
tiff("WHR_Full_model_QQ_plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## put in which phenotype you are looking at  
trait <- "Waist-to-hip ratio" #example 

## qq plot 
qq(whr_ewas$p, main = paste("Q-Q plot of", trait, "OSCA EWAS"), pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis

text(norm_x, norm_y, "\u03BB = 3.7")

## save out 
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 
# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-whr_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#[1] 3.746786


## Step 3 - make Manhattan plot 
ppi = 300
tiff("WHR_Full_model_Manhat_Plot.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  whr_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D","#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = paste("Waist-to-hip ratio (WHR)")
)

dev.off()



## Read in EWAS File 

body_fat = fread("body_fat_EWAS_full_model_filtered_CpGs.csv") 
body_fat_ewas = as.data.frame(body_fat)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 

## start saving out a tiff image for the qq plot 
tiff("Body_Fat_QQ_Plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## put in which phenotype you are looking at  
trait <- "Body fat (%)" #example 

## qq plot 
qq(body_fat_ewas$p, main = paste("Q-Q plot of", trait, "OSCA EWAS"), pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis



text(norm_x, norm_y, "\u03BB = 5")

## save out  
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 

# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-body_fat_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#[1] 4.955852


## Step 3 - make Manhattan plot 
ppi = 300
tiff("Body_Fat_Manhat.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  body_fat_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D", "#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = paste("Body fat (%)")
)


dev.off()



## Read in EWAS File 

Glucose = fread("Glucose_EWAS_full_model_filtered_CpGs.csv") 
Glucose_ewas = as.data.frame(Glucose)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 

## start saving out a tiff image for the qq plot 
tiff("Glucose_QQ_plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## put in which phenotype you are looking at  
trait <- "Glucose (mmol/L)" #example 

## qq plot 
qq(Glucose_ewas$p, main = paste("Q-Q plot of", trait, "OSCA EWAS"), pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis

text(norm_x, norm_y, "\u03BB = 1.8")

## save out
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 
# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-Glucose_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#[1] 1.834439


## Step 3 - make Manhattan plot 
ppi = 300
tiff("Glucose_Manhat_plot.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  Glucose_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D", "#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = paste("Glucose (mmol/L)")
)


dev.off()



## Read in EWAS File 

HDL_cholesterol = fread("HDL_cholesterol_EWAS_full_model_filtered_CpGs.csv") 
HDL_cholesterol_ewas = as.data.frame(HDL_cholesterol)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 
## start saving out a tiff image for the qq plot 
tiff("HDL_cholesterol_QQ_plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## put in which phenotype you are looking at  
trait <- "HDL cholesterol (mmol/L)" #example 

## qq plot 
qq(HDL_cholesterol_ewas$p, main = paste("Q-Q plot of", trait, "\nOSCA EWAS"), pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis



text(norm_x, norm_y, "\u03BB = 7.4")

## save out 
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 
# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-HDL_cholesterol_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#[1] 7.381642


## Step 3 - make Manhattan plot 
ppi = 300
tiff("HDL_cholesterol_Manhat_Plot.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  HDL_cholesterol_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D", "#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = paste("HDL cholesterol (mmol/L)")
)

dev.off()


## Read in EWAS File 

Total_cholesterol = fread("Total_cholesterol_EWAS_full_model_filtered_CpGs.csv") 
Total_cholesterol_ewas = as.data.frame(Total_cholesterol)


## Step 1 - Make QQ plot (checks for inflation of test statistics) 
## start saving out a tiff image for the qq plot 
tiff("Total_cholesterol_QQ_plot.tiff", width= 5*ppi, height=5*ppi, res=ppi)

## put in which phenotype you are looking at  
trait <- "Total cholesterol (mmol/L)" #example 

## qq plot 
qq(Total_cholesterol_ewas$p, main = paste("Q-Q plot of", trait, "\nOSCA EWAS"), pch = 18, col = "blue4", cex = 1.5, las = 1)

usr <- par("usr")

# Normalize coordinates to plot region (0 to 1)
norm_x <- usr[1] + (usr[2] - usr[1]) * 0.5 # Center of x-axis
norm_y <- usr[3] + (usr[4] - usr[3]) * 0.8 # Near top of y-axis


text(norm_x, norm_y, "\u03BB = 2.3")

## save out 
dev.off()


## Step 2 - calculate lambda (inflation factor - quantifies amount of inflation) 
# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-Total_cholesterol_ewas$p,1)

## Calculate lambda gc (λgc)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#[1] 2.32692


## Step 3 - make Manhattan plot 
ppi = 300
tiff("Total_cholesterol_Manhat_plot.tiff", width=6.75*ppi, height=6*ppi, res=ppi)

manhattan(
  Total_cholesterol_ewas,
  chr = "Chr",
  bp = "bp",
  p = "p",
  snp="Probe",
  col = c("#F8766D", "#619CFF"),
  chrlabs = NULL,
  suggestiveline = -log10(3.6e-08),
  genomewideline = -log10(3.6e-08),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  cex.axis = 1.2,
  cex.lab = 1.2,
  main = paste("Total cholesterol (mmol/L)")
)

dev.off()
