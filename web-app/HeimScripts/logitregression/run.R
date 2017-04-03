library(reshape2)

#main <- function(transformation = "raw", selectedPatientIDs = integer() ){
main <- function(transformationx = "raw", transformationy = "raw", selectedPatientIDs = integer() ){

    df1 = loaded_variables$datapoints_n0_s1 
    df2 = loaded_variables$datapoints_n1_s1

    # fetching parameters 
    if (nrow(df1) == 0) {
        stop(paste("Variable '", fetch_params$ontologyTerms$datapoints_n0$name, "' has no patients for subset 1"), sep="")
    }
    if (nrow(df2) == 0) {
        stop(paste("Variable '", fetch_params$ontologyTerms$datapoints_n1$name, "' has no patients for subset 1"), sep="")
    }

    if(is.null(transformationx)){
        transformationx = "raw"
    }

    if(is.null(transformationy)){
        transformationy = "raw"
    }

    # data preparation
    df = merge(df1, df2, by="Row.Label", sort=TRUE)
    colnames(df) <- c("patientID", "x", "y")
    xvar = colnames(df)[2]
    yvar = colnames(df)[3]


    # optional transformation
    if(transformationx == "raw"){
        # keep data unaltered 
    } else if(transformationx == "log2"){
        df$x = log2(df$x)  
    } else if(transformationx == "log10"){
        df$x = log10(df$x)  
    }


    # remove NA and infinity
    df = df[!is.infinite(df$x), ]
    df = df[!is.infinite(df$y), ]
    df = df[complete.cases(df), ]
    raw = df

    if(nrow(df)<=0 && ncol(df)<=0){
        stop("No values available for regression")
    }

    forceNormalize = length( which( df$y < 0 )) || length( which( df$y > 1 ))

    # transformation to keep y data in range [0, 1]
    if(transformationy == "raw"){
        #keep values unaltered
    }
    if(transformationy == "clamp"){
        df$y = y=min(max(df$y, 0), 1)
    }
    if((transformationy == "normalize") || ((transformationy == "raw") && (forceNormalize))){
        miny = min(df$y)
        maxy = max(df$y)

        normy = maxy - miny
        if(normy == 0){
            normy = maxy
        }
        df$y = ((df$y-miny)/normy)
    }  


    # get selected patients
    if(length(selectedPatientIDs) > 0){
        df <- df[df$patientID %in% selectedPatientIDs, ]
    }


    # regression
    inform = as.formula(paste(yvar, xvar, sep="~"))
    glm.out = glm(formula=inform, family=binomial(logit), data=df)

    # correlation of input and fitted points
    correlation = cor(glm.out$fitted.values, df$y)
    pvalue = anova(glm.out, test="Chisq")

    fitted = data.frame(x=df$x, y=glm.out$fitted.values)

    output = list(
        raw = raw,
        data = df,
        fitted = fitted,
        deviance = glm.out$deviance,
        residuals = glm.out$residuals,
        xArrLabel = fetch_params$ontologyTerms$datapoints_n0$fullName,
        yArrLabel = fetch_params$ontologyTerms$datapoints_n1$fullName,
        transformationx = transformationx,
        transformationy = transformationy,
        selected = selectedPatientIDs,
        forceNormalize = forceNormalize,
        pvalue = pvalue
    ) 
    toJSON(output)
}