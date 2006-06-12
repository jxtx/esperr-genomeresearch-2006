score = scan( commandArgs()[2] )
cpg = scan( commandArgs()[3] )
cons = scan( commandArgs()[4] )

outfile = commandArgs()[5]

tbl = as.data.frame( cbind( cpg, cons, score ) )

fit = lm( score ~ cpg + cons, data=tbl )

F = vector( "numeric", nrow( tbl ) )
for ( i in 1:nrow(tbl) ) { F[i] = fit$residuals[as.character(i)] }

write( F, outfile, ncolumns=1 )

