
suppressPackageStartupMessages(library(RPostgreSQL))
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(gsubfn))
#setwd('/Users/scottalexandermalec/Projects/semiramis')
#########
# DBsetup
#print("loading DB setup")
#########
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "causalehr",
    host = "localhost",
    port = 5432,
    user = "scottalexandermalec",
    password = "mandarin"
  )


options(width = 250)
#########

cuiLookup <- function(cuis) {
  covnames <- c()
  for (cui in cuis) {
    meaning <- dbGetQuery(con, paste("select distinct lower(covname) as covname from causalconcepts cc where cc.cui = '", cui, "'; ", sep = ""))
    #print(meaning[,1])
    covnames <- c(covnames, meaning[,1])
  }
  sort(unlist(covnames))
}
#cuiLookup(confounders)

exposure = " ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517') "
outcome = " ('C0002395', 'C0494463', 'C0276496') "

#########
# manually defined exclusions
#########
manuallyDefinedExclusions <- c('C0450442', 'C0243077', 'C0021083', 'C0013216', 'C0009429', 'C0199176', 'C0085104',
                               'C0017296', 'C0597357', 'C1519033', 'C1280903', 'C0003289', 'C0279025', 'C0279025',
                               'C0040808', 'C0009244', 'C0175723', 'C0011581', 'C0026868', 'C0683465', 'C1519595',
                               'C0042890', 'C0150321', 'C0581602', 'C1874451', 'C0012854', 'C0033972', 'C0279493',
                               'C0014935', 'C0086132', 'C0344315', 'C0016452', 'C1137094', 'C2718059', 'C0013806',
                               'C0183185', 'C0678908', 'C2827774', 'C0031765', 'C0003827', 'C0034991', 'C0581601',
                               'C0282402', 'C0150593', 'C1254359', 'C3858690', 'C1096024', 'C0020431', 'C1293131',
                               'C1964029', 'C0204727', 'C0005893', 'C0185027', 'C0439775', 'C0418981', 'C0152060',
                               'C0183210', 'C0870230', 'C0677850', 'C0376626', 'C0677850', 'C0949766', 'C0150133',
                               'C1515119', 'C0079613', 'C0455142', 'C0585941', 'C0279494', 'C0920425', 'C0033968', 
                               'C0935576', 'C1527374', 'C1269683', 'C0034618', 'C0376547', 'C0020431', 'C0185117',
                               'C0204514', 'C0439857', 'C0680038', 'C0846672', 'C0015618', 'C0182537', 'C0183115', 
                               'C0220908', 'C0279494', 'C0935576', 'C1268930', 'C0013218', 'C0034764', 'C0183210', 
                               'C0232164', 'C0281481', 'C1139730', 'C1511300', 'C1512629', 'C2986605', 'C3177188',
                               'C0030567', 'C0497327', 'C0474169', 'C0011269', 'C0599917', 'C0599917')

length(manuallyDefinedExclusions)

#########
# Covariates that are measured for subjects in sample from electronic health database
#########
getGoodCovars <- function(n) {
  fNames <- list.files(path = "data/r3dataset/data")
  i = 0
  goodData <- c()
  for (f in fNames) {
    if (sum(read.csv(file = paste("data/r3dataset/data/", f, sep = ""), header = FALSE)) >= n) 
    { #print(f) 
      i = i + 1 
      #print(i)
      goodData <- c(goodData, f) 
    }
  }
  #print(goodData)
  goodData <- gsubfn(goodData, pattern = "\\.csv", replacement = "")
  return(goodData)
}

#predicates <- " ('CAUSES', 'PREVENTS', 'PREDISPOSES') " 
# NOTE: 'TREATS' -- excluding because of the ambiguous causal meaning of the verb 'to treat'
predicates <- " ('CAUSES', 'PREDISPOSES') " #, 'PREVENTS', 'AFFECTS', 'DISRUPTS', 'AUGMENTS', 'STIMULATES', 'INHIBITS') " # 'AFFECTS', 'DISRUPTS', 'AUGMENTS', 'TREATS',
#predicates <- " ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS', 'AFFECTS', 'DISRUPTS', 'AUGMENTS', 'STIMULATES', 'INHIBITS') "

#########
# get confounders
#########
getConfounders <- function(predicates) {
  confounders <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar, COUNT(*) AS scnt, cp.subject_name --, string_agg(cp.subject_semtype) as subject_semtype
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ", exposure, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.subject_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.subject_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.subject_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY subject_cui, subject_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name --, string_agg(cp.subject_semtype) as subject_semtype
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ", outcome, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.subject_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.subject_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.subject_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, SUM(zx.scnt) as theta1, SUM(zy.scnt) AS theta2, replace(lower(zx.subject_name), ' ', '_') as covarname --, string_agg(zx.subject_semtype) 
                                       FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar 
                                       GROUP BY ZX.covar, replace(lower(zx.subject_name), ' ', '_')
                                       ORDER BY theta2 desc, theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawConfounders_for_depression2AD.txt", x = confounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  confounders <- subset(confounders, confounders[,c(2)] > 2)
  confounders <- na.omit(subset(confounders, confounders[,c(3)] > 2))
  removeList <- c('uptake', 'disability_nos', 'cell_physiology', 'signal_transduction_pathways', 'emotional', 'vision', 'growth', 'energy_metabolism', 'dependence', 'signs_and_symptoms', 'growth_factor', 'indicators', 'equilibrium', 'sedentary', 'drugs,_investigational', 'food', 'genetic_aspects', 'basis', 'oxygen', 'procedure_findings:finding:point_in_time:^patient:narrative', 'biological_response_modifiers', 'gene_polymorphism', 'transcript', 'chromosome_pairing', 'carrier_of_disorder', 'metric','derivatives', 'small_molecule', 'ingredient', 'rna,_circular', 'atrophic', 'antagonists', 'agent', 'memantine', 'disease_progression', 'micrornas', 'imbalance', 'study_models', 'receptor', 'amyloid', 'app_gene', 'hypersensitivity', 'possible', 'risk_factors', 'substance', 'nutrients', 'impact_gene', 'uncertainty', 'personality_traits', 'instability', 'sedentary', 'fruit', 'serotonin', 'discrepancy')
  removeList2 <- c('rna,_untranslated', 'autophagy', 'predicted', 'vitamins', 'app_gene', 'memantine', 'amyloid', 'galantamine', 'alzheimer\'s_disease','therapeutic_agent_(substance)', 'glutamate_receptor', 'cannabinoid_receptor')
  confounderList <- setdiff(confounders[,4], removeList)
  confounderList <- setdiff(confounderList, removeList2)
  confounders <- na.omit(subset(confounders, confounders[,4] %in% confounderList))
  write.table(file = "covariates/rawConfounders_for_depression2AD.txt", x = confounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(confounders)
}


confounders <- getConfounders(predicates)
print(confounders)
sort(confounders[,4])
#########
#
#########
#

getColliders <- function(predicates) { ########## FIX SEMTYPE
  colliders <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ", exposure, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.object_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.object_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.object_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY object_cui, object_name),
                        ZY AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ", outcome, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.object_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.object_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.object_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY cp.object_cui, cp.object_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, replace(lower(zx.object_name), ' ', '_') as covarname FROM ZX AS zx, ZY AS zy 
                                     WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  colliders <- subset(colliders, colliders[,c(2)] > 1)
  colliders <- subset(colliders, colliders[,c(3)] > 1)
  removeList <- c('uptake', 'disability_nos', 'cell_physiology', 'signal_transduction_pathways', 'emotional', 'vision', 'growth', 'energy_metabolism', 'dependence', 'signs_and_symptoms', 'growth_factor', 'indicators', 'equilibrium', 'sedentary', 'drugs,_investigational', 'food', 'genetic_aspects', 'basis', 'oxygen', 'procedure_findings:finding:point_in_time:^patient:narrative', 'biological_response_modifiers', 'gene_polymorphism', 'transcript', 'chromosome_pairing', 'carrier_of_disorder', 'metric','derivatives', 'small_molecule', 'ingredient', 'rna,_circular', 'atrophic', 'antagonists', 'agent', 'memantine', 'disease_progression', 'micrornas', 'imbalance', 'study_models', 'receptor', 'amyloid', 'app_gene', 'hypersensitivity', 'possible', 'risk_factors', 'substance', 'nutrients', 'impact_gene', 'uncertainty', 'personality_traits', 'instability', 'sedentary', 'fruit', 'serotonin', 'discrepancy')
  removeList2 <- c('rna,_untranslated', 'autophagy', 'predicted', 'vitamins', 'app_gene', 'memantine', 'amyloid', 'galantamine','alzheimer\'s_disease', 'therapeutic_agent_(substance)', 'glutamate_receptor', 'cannabinoid_receptor')
  colliderList <- setdiff(colliders[,4], removeList)
  colliderList <- setdiff(colliderList, removeList2)
  colliders <- na.omit(subset(colliders, colliders[,4] %in% colliderList))
  write.table(file = "covariates/rawColliders_for_depression2AD.txt", x = colliders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(colliders)
}

#getColliders(predicates)
 
colliders <- getColliders(predicates)
colliders
sort(colliders[,4])

getMediators <- function(predicates) {
  mediators <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name as covarname
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ", exposure, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.object_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.object_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.object_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY object_cui, object_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ", outcome, "
                        AND cp.predicate IN ", predicates, "
                        AND cp.subject_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.subject_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.subject_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.covarname) as covarname FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawMediators_for_depression2AD.txt", x = mediators, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  mediators <- subset(mediators, mediators[,c(2)] > 1)
  mediators <- subset(mediators, mediators[,c(3)] > 1)
  removeList <- c('uptake', 'disability_nos', 'cell_physiology', 'signal_transduction_pathways', 'emotional', 'vision', 'growth', 'energy_metabolism', 'dependence', 'signs_and_symptoms', 'growth_factor', 'indicators', 'equilibrium', 'sedentary', 'drugs,_investigational', 'food', 'genetic_aspects', 'basis', 'oxygen', 'procedure_findings:finding:point_in_time:^patient:narrative', 'biological_response_modifiers', 'gene_polymorphism', 'transcript', 'chromosome_pairing', 'carrier_of_disorder', 'metric','derivatives', 'small_molecule', 'ingredient', 'rna,_circular', 'atrophic', 'antagonists', 'agent', 'memantine', 'disease_progression', 'micrornas', 'imbalance', 'study_models', 'receptor', 'amyloid', 'app_gene', 'hypersensitivity', 'possible', 'risk_factors', 'substance', 'nutrients', 'impact_gene', 'uncertainty', 'personality_traits', 'instability', 'sedentary', 'fruit', 'serotonin', 'discrepancy')
  removeList2 <- c('rna,_untranslated', 'autophagy', 'predicted', 'vitamins', 'app_gene', 'memantine', 'amyloid', 'galantamine', 'alzheimer\'s_disease', 'therapeutic_agent_(substance)', 'glutamate_receptor', 'cannabinoid_receptor')
  mediatorList <- setdiff(mediators[,4], removeList)
  mediatorList <- setdiff(mediatorList, removeList2)
  mediators <- na.omit(subset(mediators, mediators[,4] %in% mediatorList))
  write.table(file = "covariates/rawMediatorss_for_depression2AD.txt", x = mediators, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(mediators)
}

mediators <- getMediators(predicates)
sort(mediators[,4])

# getPrecisionVars <- function(predicates) {
#   precisionVars <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
#                         FROM causalpredications cp, citations c
#                         WHERE cp.pmid = c.pmid
#                         AND CAST(pyear AS NUMERIC) >= 2010
#                         AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
#                         AND cp.predicate IN ", predicates, "
#                         GROUP BY cp.subject_cui, cp.subject_name)
#                         SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
#   write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
#   return(precisionVars)
# }
# 
# precisionVars <- getPrecisionVars(predicates)

# getOutcomeEffects <- function(predicates) {
#   outcomeEffects <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name as covarname
#                         FROM causalpredications cp, citations c
#                         WHERE cp.pmid = c.pmid
#                         AND CAST(pyear AS NUMERIC) >= 2010
#                         AND cp.subject_cui IN ('C0002395', 'C0494463', 'C0276496')
#                         AND cp.predicate IN ", predicates, "
#                         GROUP BY cp.object_cui, cp.object_name)
#                         SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
#   #write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
#   return(outcomeEffects)
# }
# 
# outcomeEffects <- getOutcomeEffects(predicates)


# getExposureEffects <- function(predicates) {
#   exposureEffects <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name as covarname
#                         FROM causalpredications cp, citations c
#                         WHERE cp.pmid = c.pmid
#                         AND CAST(pyear AS NUMERIC) >= 2010
#                         AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
#                         AND cp.predicate IN ", predicates, "
#                         GROUP BY cp.object_cui, cp.object_name)
#                         SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
#   #write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
#   return(exposureEffects)
# }
# exposureEffects <- getExposureEffects(predicates)

#confounders[,4]
#colliders[,4]
#mediators[,4]
#precisionVars[,3]

confounderList <- setdiff(setdiff(confounders[,4], colliders[,4]), mediators[,4])

confoundersAlmostCooked <- subset(confounders, confounders[,4] %in% confounderList)

confs <- unlist(confoundersAlmostCooked[,1])
#conf1 <- confs[1]
#conf2 <- confs[2]
#print(c1)
#print(c2)

getPotentialMBiasVars <- function(predicates, c1, c2) { ########## FIX SEMTYPE
  #if (file.exists(!paste("covariates/mBiasVarsMDD2AD-", c1, "-", c2, ".txt", sep = ""))) {
  mbiasSQL <- paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name, cp.subject_name as c1
                        FROM causalpredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('", c1, "')
                        AND cp.predicate IN ", predicates, "
                        AND cp.object_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.object_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.object_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY cp.object_cui, cp.object_name, cp.subject_name),
                        ZY AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name, cp.subject_name as c2
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('", c2, "')
                        AND cp.predicate IN ", predicates, "
                        AND cp.object_semtype NOT IN ('orgm','bpoc','diap','ortf','bsoj','inpo','tisu','topp','mamm','inpr','geoa','hlca','bdsy','blor','hcro','lbpr','inbe','orga','menp','mnob','humn','amph','plnt','spco','anim','resa','anab','eehu','tmco','edac','ftcn','ocdi', 'dora','qnco','orgt','npop','qlco','podg','prog','bird','mcha','rept','fish','phob','socb','idcn','popg', 'bmod','emod','aggp','famg','rnlw','mbrt','pros','lang','ocac','gora','medd')
                        AND cp.object_cui NOT IN ('C0001687','C0002526','C0003043','C0003062','C0005515','C0009566','C0012634','C0013227','C0021521','C0021948','C0027361','C0027362','C0027363','C0028622','C0029224','C0029235','C0030705','C0039082','C0039796','C0087111','C0159028','C0178310','C0178341','C0178353','C0178355','C0178359','C0243192','C0422820','C0424450','C0436606','C0442826','C0476466','C0478681','C0478682','C0480773','C0481349','C0481370','C0557587','C0565657','C0580210','C0580211','C0589603','C0596048','C0596090','C0597010','C0597237','C0597240','C0677042','C0687732','C1257890','C1258127','C1318101','C1457887','C1609432','C0007634','C0020114','C0237401','C0011900','C1273869','C0449851','C0277785','C0184661','C1273870','C0185125','C0879626','C0004927','C0936012','C0311392','C0597198','C0018684','C0042567','C0029921','C0683971','C0016163','C0024660','C0687133','C0037080','C0680022','C1185740','C0871261','C0544461','C1260954','C0877248','C0242485','C0205147','C0486805','C0005839','C0021562','C0205148','C0031843','C0040223','C0205145','C0205400','C0086388','C0014406','C0520510','C0035168','C0029237','C0277784','C0001779','C0542559','C0035647','C0025664','C0700287','C0678587','C0205099','C0205146','C0237753','C0441800','C0449719','C0348026','C0008902','C0586173','C0332479','C0807955','C0559546','C0031845','C0678594','C0439792','C0557854','C1522240','C1527144','C0449234','C0542341','C0079809','C0205094','C0037455','C0025118','C0441471','C0441987','C0439534','C0392360','C0456603','C0699733','C0036397','C0725066','C0496675','C0282354','C0015127','C1273937','C1368999','C0442804','C0449286','C0205082','C0814472',    'C1551338','C0599883','C0450429','C1299582','C0336791','C0443177','C0025080','C1372798','C0028811','C0205246','C0449445','C0332185','C0332307','C0443228','C1516635','C0376636','C0221423','C0037778','C0199168','C0008949','C0014442'    ,'C0456387','C1265611','C0243113','C0549177','C0229962','C0600686','C1254351','C0243095','C1444647','C0033684','C0338067','C0441712','C0679607','C0808233','C1373236','C0243082','C1306673','C1524062','C0002085','C0243071','C0238767','C0005508','C0392747','C0008633','C0205195','C0205198','C0456205','C0521116','C0011155','C1527240','C1527148','C0743223','C0178602','C1446466','C0013879','C0015295','C1521761','C1522492','C0017337','C0017428','C0017431','C0079411','C0018591','C0019932','C0021149','C0233077','C0021920','C0022173','C1517945','C0680220','C0870883','C0567416','C0596988','C0243132','C0029016','C1550456','C0243123','C0030956','C0851347','C0031328','C0031327','C00314','C1514468','C0033268','C0449258','C0871161','C1521828','C0443286','C1547039','C1514873','C0035668','C0439793','C0205171','C0449438','C1547045','C0449913','C0042153','C0205419','C1441526','C1140999','C0679670','C0431085','C1185625','C1552130','C1553702','C1547020','C0242114','C0439165','C0679646','C0599755','C0681850','6275','6285')
                        AND cp.object_name NOT IN ('inhibitors', 'Single Nucleotide Polymorphism', 'multiple pathologies', 'Natural Products', 'TRANSCRIPTION FACTOR', 'Transcription, Genetic', 'Molecular Target', 'Gene Expression', 'Signal Transduction', 'Homeostasis', 'Accumulation', 'Syndrome', 'Disease', 'Genetic disorders', 'Polymorphism, Genetic', 'Mental disorders', 'Dementia', 'DNA', 'Alzheimer disease, familial, type 3', 'Pathogenesis')
                        GROUP BY cp.object_cui, cp.object_name, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, replace(lower(zx.object_name), ' ', '_') as covarname, replace(lower(zx.c1), ' ', '_') as c1, replace(lower(zy.c2), ' ', '_') as c2 FROM ZX AS zx, ZY AS zy 
                                     WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = "")
  mbiasVars <- dbGetQuery(con, mbiasSQL)
  mbiasVars <- subset(mbiasVars, mbiasVars[,c(2)] > 1)
  mbiasVars <- subset(mbiasVars, mbiasVars[,c(3)] > 1)
  removeList <- c('uptake', 'disability_nos', 'cell_physiology', 'signal_transduction_pathways', 'emotional', 'vision', 'growth', 'energy_metabolism', 'dependence', 'signs_and_symptoms', 'growth_factor', 'indicators', 'equilibrium', 'sedentary', 'drugs,_investigational', 'food', 'genetic_aspects', 'basis', 'oxygen', 'procedure_findings:finding:point_in_time:^patient:narrative', 'biological_response_modifiers', 'gene_polymorphism', 'transcript', 'chromosome_pairing', 'carrier_of_disorder', 'metric','derivatives', 'small_molecule', 'ingredient', 'rna,_circular', 'atrophic', 'antagonists', 'agent', 'memantine', 'disease_progression', 'micrornas', 'imbalance', 'study_models', 'receptor', 'amyloid', 'app_gene', 'hypersensitivity', 'possible', 'risk_factors', 'substance', 'nutrients', 'impact_gene', 'uncertainty', 'personality_traits', 'instability', 'sedentary', 'fruit', 'serotonin', 'discrepancy')
  removeList2 <- c('rna,_untranslated', 'autophagy', 'predicted', 'vitamins', 'app_gene', 'memantine', 'amyloid', 'alzheimer\'s_disease', 'galantamine', 'therapeutic_agent_(substance)', 'glutamate_receptor', 'cannabinoid_receptor')
  mbiasVarList <- setdiff(mbiasVars[,4], removeList)
  mbiasVarList <- setdiff(mbiasVarList, removeList2)
  mbiasVars <- na.omit(subset(mbiasVars, mbiasVars[,4] %in% mbiasVarList))
  return(mbiasVars)
}

#mbiasVars <- getPotentialMBiasVars(predicates, confs[1], confs[2])
#getPotentialMBiasVars(predicates, "C0243077", "C0033968")

confoundersNotMbias <- confoundersAlmostCooked[,4]

for (conf1 in confs) {
  for (conf2 in confs) {
    if ((conf1 != conf2) && (!file.exists(paste("covariates/mBiasVarsMDD2AD-", conf2, "-", conf1, ".txt", sep = "")))  && (!file.exists(paste("covariates/mBiasVarsMDD2AD-", conf1, "-", conf2, ".txt", sep = "")))) {
      print(paste("RUNNING mbias variable search - conf1: ", conf1, " conf2: ", conf2, sep = ""))
      mbiasVars <- getPotentialMBiasVars(predicates, conf1, conf2)
      print(mbiasVars)
      print("###############")
      print("BINGOS")
      mbiasVariables <- intersect(mbiasVars[,4], confoundersAlmostCooked[,4]) # list
      confoundersNotMbias <- setdiff(confoundersNotMbias, mbiasVariables)
      print("###############")
      print("############### NOT M-BIAS Variables (so far)")
      print(confoundersNotMbias)
      cleanConfounders <- subset(confoundersAlmostCooked, confoundersAlmostCooked[,4] %in% confoundersNotMbias) # new
      print(nrow(cleanConfounders)) # new
      print("###############")
      print("###############")
      print("###############")
      print("mbiasVars")
      mbiasConfounders <- na.omit(subset(mbiasVars, mbiasVars[,4] %in% mbiasVariables)) # table
      ######
      print("###############")
      print("mbiasconfounders")
      print(mbiasConfounders)
      print("###############")
      write.table(file = paste("covariates/mBiasVarsMDD2AD-", conf1, "-", conf2, ".txt", sep = ""), x = mbiasConfounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
      #return(colliders)
    }
  }
}


confoundersCooked <- setdiff(setdiff(cleanConfounders[,4], colliders[,4]), mediators[,4])
print(confoundersNotMbias)
write.table(file = "covariates/cookedCleanConfounders_for_AD.txt", x = confoundersCooked, row.names = FALSE, col.names = FALSE, quote = FALSE)  

#print(confoundersCooked)
#confoundersCookedCooked <- setdiff(confoundersCooked, outcomeEffects[,3])
#confoundersCookedCooked

#print(colliders)
collidersCooked <- setdiff(setdiff(colliders[,4], confoundersCooked[,4]), mediators[,4])
print(collidersCooked)
write.table(file = "covariates/cookedColliders_for_AD.txt", x = collidersCooked, row.names = FALSE, col.names = FALSE, quote = FALSE)  

print(mediators)
mediatorsCooked <- setdiff(setdiff(mediators[,4], collidersCooked[,4]), confoundersCooked[,4])
print(mediatorsCooked)
write.table(file = "covariates/cookedMediators_for_AD.txt", x = mediatorsCooked, row.names = FALSE, col.names = FALSE, quote = FALSE)  

chimeras <- intersect(intersect(confoundersCooked[,4], collidersCooked[,4]), mediatorsCooked[,4])
write.table(file = "covariates/cookedChimeras_for_MDD2AD.txt", x = chimeras, row.names = FALSE, col.names = FALSE, quote = FALSE)  

