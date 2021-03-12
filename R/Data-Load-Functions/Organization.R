###################################
# Longitudinally Ordering Samples # 
###################################

# Creating a set of lists of the order of the samples in each cohort. 
TLML_1_gDNA <- c("apheresis_2013_9", 'infusion_2013_8', "FU_01_2014_1", "FU_02_2014_4", "FU_03_2014_7")
TLML_1_cDNA <- c('4W_2013_11', 'FU_07_2015_2', 'FU_10_2015_6', 'FU_12_2015_10', 'FU_13_2015_11', 'FU_16_2016_6')
TLML_4_gDNA <- c('baseline_2015_10', 'infusion_2013_11', '4W_2014_12', 'FU_01_2015_1')
TLML_4_cDNA <- c('baseline_2015_10', 'infusion_2013_11', '4W_2014_12', 'FU_01_2015_1')
TLML_4_cfDNA <- c('baseline_2015_10', '4W_2014_12')
TLML_7_gDNA <- c('apheresis_2014_3', 'infusion_2014_2', '4W_2014_4', '12W_2014_6', '17W_2014_7')
TLML_7_cDNA <- c('apheresis_2014_3', 'infusion_2014_2', '4W_2014_4', '12W_2014_6')
TLML_16gDNA <- c('baseline_2015_1', 'infusion_2014_4', '4W_2015_2', 'FU_01_2015_3', 'FU_02_2015_5', 'FU_03_2015_6', 'FU_04_2015_7', 'FU_05_2015_9')
TLML_16cDNA <- c('baseline_2015_1', 'infusion_2014_4', 'FU_01_2015_3', 'FU_02_2015_5', 'FU_03_2015_6', 'FU_04_2015_7', 'FU_05_2015_9')
TLML_16cfDNA <- c('baseline_2015_1', '4W_2015_2')
TLML_18gDNA <- c('apheresis_2014_8', 'infusion_2014_6', '4W_2014_9')
TLML_18cDNA <- c('apheresis_2014_8', 'infusion_2014_6', '4W_2014_9')
TLML_20gDNA <- c('baseline_2_2015_1', 'infusion_2014_7', '4W_2015_2', 'FU_01_2015_3', 'FU_02_2015_6')
TLML_20cDNA <- c('baseline_2_2015_1', 'infusion_2014_7', '4W_2015_2', 'FU_01_2015_3', 'FU_02_2015_6')
TLML_20cfDNA <- c('baseline_2_2015_1', '4W_2015_2')
TLML_22gDNA <- c('baseline_2014_12', 'infusion_2014_9', '4W_2015_1', 'FU_01_2015_2', 'FU_02_2015_3', 'FU_03_2015_4')
TLML_22cDNA <- c('infusion_2014_9', '4W_2015_1', 'FU_01_2015_2', 'FU_02_2015_3', 'FU_03_2015_4')
TLML_22cfDNA <- c('baseline_2014_12', '4W_2015_1')
TLML_26gDNA <- c('apheresis_2014_11', 'infusion_2014_6','4W_2015_3', 'FU_01_2015_7', 'FU_02_2015_9')
TLML_26cDNA <- c('baseline_2015_1', '4W_2015_3')
TLML_26cfDNA <- c('baseline_2015_1', '4W_2015_3')
TLML_29gDNA <- c('baseline_2016_1', 'infusion_2015_10', '4W_2016_2')
TLML_29cDNA <- c('baseline_2016_1', 'infusion_2015_10', '4W_2016_2')
TLML_29cfDNA <- c('baseline_2016_1', '4W_2016_2')
TLDC_1_gDNA <- c('baseline_2015_11', 'infusion_NA_NA', '4W_2015_12')
TLDC_1_cDNA <- c('baseline_2015_11', 'infusion_NA_NA', '4W_2015_12')
TLDC_1_cfDNA <- c('baseline_2015_11', '4W_2015_12')

#####################
# Loading Libraries #
#####################

# Loading required libraries & options (if not previously installed, you must first install these packages)
library(ggplot2)
library(ggalluvial)
library(immunarch)
library(readxl)
library(cowplot)
options(scipen = 999)

#########################################
# Loading Manual Distinct Color Pallete #
#########################################

# Maximally distant color pallete excluding black and white

colors <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6",
           "#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF",
           "#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF",
           "#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF",
           "#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1",
           "#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD",
           "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700",
           "#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700",
           "#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55",
           "#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489",
           "#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5",
           "#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1",
           "#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7",
           "#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78",
           "#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7",
           "#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966",
           "#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891",
           "#15A08A","#BC65E9","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527",
           "#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199",
           "#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03",
           "#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42",
           "#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064",
           "#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D",
           "#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25",
           "#252F99","#00CCFF","#674E60","#FC009C","#92896B")