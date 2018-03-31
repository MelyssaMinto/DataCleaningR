IsDate <- function(mydate, date.format = "%d/%m/%y") {
    # Check if field is a date using as.Date that looks for unambiguous dates
    # Input: mydate - suspected date and optional date format string
    # Output: TRUE if thinks it is a date
    tryCatch(!is.na(as.Date(mydate, date.format)),  
             error = function(err) {FALSE})  
}

insertNA <- function(ref, new, v){
  # This function will insert an NA if missing from input file. It checks if the difference between entries 
  # are the same. If not, it iserts NA. 
  # Input: ref - reference index
  #       new - index of interest
  # output v - vector of interest
  if(length(ref)==length(new)) return(v)
  
  diff1<-diff(ref)
  diff2<-diff(new)
  
  for(i in 1:length(diff1)){
    if(diff1[i]!=diff2[i]){
      diff2 = c(diff2[1:(i-1)], diff1[i], diff2[i] - diff1[i], diff2[(i+1):length(diff2)])
      v = c(v[1:(i-1)],NA, v[i:length(v)])
    }
  }
  return(v)
}

extract<-function(df) {
  #each row has one element, and each value is proceeded by the values name. 
  #So to extract the data this function will search for the values' names and then pull 
  #the following row. For example if the word smiles is found on rows 1,5,and 9 then the 
  #function will pull rows 2,6, and 10 asthe smiles ID.
  
  ###extracting smiles
  s = grep("smiles", df$V1, ignore.case = T)
  smiles = df$V1[s+1]

  ####extracting CASRN
  c = grep("casrn", df$V1, ignore.case = T)
  #some casrns have been converted to data format MM/DD/YYYY so, that will have to be reformatted 
  #by changing "/" to "-" ===>  MM/DD/YYYY to MM-DD-YYYY
  #and swithing the month and year positions ad  ==> YYYY-DD-MM 
  casrn = as.character(df$V1[c+1])
  casrn = ifelse(IsDate(casrn), format(as.Date(casrn, format ="%m/%d/%Y" ), "%Y-%d-%m"), casrn)
  casrn = insertNA(s, c, casrn)
                      
  ###extracting DSSToxid 
  d = grep("dsstox_substance_id", df$V1, ignore.case = T)
  dsstoxID = as.character(df$V1[d+1])
  dsstoxID = insertNA(c, d, dsstoxID)
  
  ###extracting Molecular Weight
  mw = grep("MW", df$V1, ignore.case = T)
  MW = as.numeric(as.character(df$V1[mw+1]))
  MW = insertNA(c, mw, MW)
  
  ###extracting batch
  b = grep("batch", df$V1, ignore.case = T)
  batch = as.numeric(as.character(df$V1[b+1]))
  batch = insertNA(c, b, batch)
  
  ###extract Cramers rules
  cr = grep("Cramer rule", df$V1, ignore.case = T)
  CramerRules = as.character(df$V1[cr+1])
  CramerRules = insertNA(c, cr, CramerRules)
  
  ###Extract Cramer tree
  ct = grep("CramerTree", df$V1, ignore.case = T)
  CramerTree = as.character(df$V1[ct+1])
  CramerTree = insertNA(c, ct, CramerTree)
  
  ###Extract Kreos Tree
  kt =  grep("<Kroes TTC decision tree>", df$V1, ignore.case = T)
  KreosTree = as.character(df$V1[kt+1])
  KreosTree = insertNA(c, kt, KreosTree)
  
  ###Extract Kreos Tree Explanation
  kte =  grep("<Kroes TTC decision tree#", df$V1, ignore.case = T)
  KreosTreeExp = as.character(df$V1[kte+1])
  KreosTreeExp = insertNA(c, kte, KreosTreeExp)
  
  data<-data.frame(smiles, 
                   casrn, 
                   dsstoxID, 
                   MW, 
                   batch, 
                   CramerRules, 
                   CramerTree, 
                   KreosTree, 
                   KreosTreeExp,
                   stringsAsFactors = F)
  
  return(data)
}