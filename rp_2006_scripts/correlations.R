tbl = read.delim( "F2.dat", header=T )

options( width=150 ) 

classes = unique( sort( as.character( tbl$dataset ) ) )

cols = c( "phastCons", "seq_cpg", "seq_gc", "rp", "prom", "distal", "F" )

cat( "\n---- Pairwise complete ------------------------------------------------\n" )

for ( c in classes )
{
    cat( "\nDataset ", c, "\n\n" )

    t = subset( tbl, dataset == c, select=cols )
    print( cor( t, use="pairwise.complete.obs" ) )
}

cat( "\n---- All complete -----------------------------------------------------\n" )

for ( c in classes )
{
    cat( "\nDataset ", c, "\n\n" )

    t = subset( tbl, dataset == c, select=cols )
    print( cor( t, use="complete.obs" ) )
}
