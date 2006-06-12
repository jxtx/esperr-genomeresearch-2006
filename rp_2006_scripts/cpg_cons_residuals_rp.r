score = scan( commandArgs()[2] )
gc = scan( commandArgs()[3] )
cons = scan( commandArgs()[4] )

outfile = commandArgs()[5]

# tbl = as.data.frame( cbind( cpg, cons, score ) )
# fit = lm( score ~ cpg + cons, data=tbl )
# F = vector( "numeric", nrow( tbl ) )
#for ( i in 1:nrow(tbl) ) { F[i] = fit$residuals[as.character(i)] }

## Using aligned CpG content
# F = score - ( 0.04391 + 1.12675 * cpg + 0.04695 * cons )

## Using GC content
# F = score - ( -0.18792 + 0.05637 * cons + 0.46588 * gc )

# Using GC content on REG + AR
F = score - ( -0.34866 + 0.17174 * cons + 0.63263 * gc )

## Using GC content on encode NCNR
#F = score - ( -0.268657 + 0.475585 * gc + 0.144657 * cons )

write( F, outfile, ncolumns=1 )
