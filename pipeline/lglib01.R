
args = commandArgs(trailingOnly = TRUE)
options (width=300)
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, n=10,m=10) {
	if (is.null (dim (data))) {
		size = length
		if (length (data) < 10) n = length(data)
	}else {
		size = dim
		if (nrow (data) < 10) n = nrow(data)
		if (ncol (data) < 10) m = ncol(data)
	}

	filename = deparse (substitute (data))
	message (paste (deparse (substitute (data)),":  "))
	print (size (data))
	if (is.null (dim (data)))
		print (data [1:m])
	else if (ncol (data) < 10) 
		print (data[1:n,])
	else if (nrow (data) < 10)
		print (data[,1:m])
	else 
		print (data [1:n, 1:m])

	write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
	
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
