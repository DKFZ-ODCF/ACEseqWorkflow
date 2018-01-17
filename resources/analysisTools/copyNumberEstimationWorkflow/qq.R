# == title
# Simple template for text interpolation
#
# == param
# -text         template text
# -replacemmet  a list containing variable name and values. If it is not specified,
#    variables in parent frame will be used.
# -code.pattern patterns of the code markings. Default is @{CODE}, which means
#    CODE part will be replaced with its real value
#
# == details
# I like text interpolation in Perl. But in R, if you want to connect plain text and variables,
# you need to use `paste` or functions in `stringr`. However, if there are so many variables or 
# or many quotes in the string you want to construct, it would be a little terrible.
#   
# So, this function allows you to construct strings as Perl style.
#
# Name for this function means double quote.
qq = function(text, replacement = parent.frame(), code.pattern = NULL) {
 
    if (is.null(code.pattern)) {
        code.pattern = "\\@\\{CODE\\}"
    }
    if(length(text) != 1) {
        stop("length of the text should be 1.")
    }
 
    lines = strsplit(text, "\n")[[1]]
    if(length(lines) == 0) {
        lines = ""
    }
    newlines = character(length(lines))
 
    # import variables in replacement
    # bug: need to new an environment and assign the value in this envrionment
    # and specify the environment in 'eval'
    
	if(!is.null(replacement)) {
		if(is.environment(replacement)) {
			e = replacement
		} else {
			e = new.env()
			for(varname in names(replacement)) {
				assign(varname, replacement[[varname]], envir = e)   
			}
		}
    } else {
        e = .GlobalEnv
    }

    for (i in 1:length(lines)) {
 
        # check wether there are code replacements
        code = find_code(code.pattern, lines[i])
        code.template = code[[1]]
        code.variable = code[[2]]
 
        if(length(code.template)) {
 
            # if there is code replacement
            # replace the code with its value
            code.result = lapply(code.variable, function(code) eval(parse(text = code), envir = e))  # anony function is the first level parent
 
            # length of the return value
            v.lines = sapply(code.result, function(x) length(x))
 
            if(max(v.lines) > 1) {
                current.line = rep(lines[i], max(v.lines))
                for(ai in 1:max(v.lines)) {
                    for(iv in 1:length(code.template)) {
                        current.line[ai] = gsub(code.template[iv],
                        code.result[[iv]][(ai-1) %% length(code.result[[iv]]) + 1],
                        current.line[ai], fixed = TRUE)
                    }
                }
                newlines[i] = paste(current.line, collapse = "\n")
            }
            else if(max(v.lines == 1)) {
                current.line = lines[i]
                for(iv in 1:length(code.template)) {
                    current.line = gsub(code.template[iv], code.result[[iv]],
                                   current.line, fixed = TRUE)
                }
                newlines[i] = current.line
            }
            else {
                newlines[i] = ""
            }
        }
        else {
            newlines[i] = lines[i]
        }
    }
 
    return(paste(newlines, collapse="\n"))
}
 
find_code = function(m, text) {
 
    if(length(text) != 1) {
        stop("text must be length of 1.")
    }
 
    m2 = gsub("CODE", ".+?", m)
 
    reg = gregexpr(m2, text, perl = TRUE)[[1]]
    v1 = character(0)
    if(reg[1] > -1) {
        v1 = sapply(1:length(reg), function(i)substr(text, as.numeric(reg)[i], as.numeric(reg)[i]+ attr(reg, "match.length")[i] - 1))
    }
    edge = strsplit(m, "CODE")[[1]]
    v2 = gsub(paste("^", edge[1], "|", edge[2], "$", sep=""), "", v1)
    return(list(template=v1, variable=v2))
}


qqcat = function(..., replacement = parent.frame()) {
	cat(qq(..., replacement = replacement), "\n")
}
