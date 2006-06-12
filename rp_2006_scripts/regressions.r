library(stepfun)

dataset = commandArgs()[2]

tbl = read.delim( dataset, comment.char="", header=T )

tbl = tbl[,4:ncol(tbl)]

cat( "\n----- Correlations ------------------------------------------------\n\n" )

print( cor( tbl, use="pairwise.complete.obs" ) )

cat( "\n----- RP regression -----------------------------------------------\n" )

print( summary( lm( rp ~ phastCons + cpg + distal, data=tbl ) ) )

cat( "----- Promoter score regression -----------------------------------\n" )

print( summary( lm( prom ~ phastCons + cpg + rp + distal, data=tbl ) ) )

cat( "----- Distal score regression -------------------------------------\n" )

print( summary( lm( distal ~ phastCons + cpg + rp + prom, data=tbl ) ) )
