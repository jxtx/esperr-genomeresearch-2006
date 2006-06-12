tbl = read.delim( "F2.dat", header=T )
# Filter out rows with NA
t = tbl[ (!is.na(tbl$seq_gc)) & (!is.na(tbl$phastCons)) & (!is.na(tbl$rp)) & (tbl$dataset != "reg_tile"), ]

# Class labels
c = vector( "numeric", nrow( t ) )
c[ t$dataset == "pos_enc_proms_tile" & t$support <5  ] = 1
c[ t$dataset == "pos_enc_proms_tile" & t$support >= 5 ] = 2
c[ t$dataset == "pos_enc_proms_tile" & t$support >= 15 ] = 3
c[ t$dataset == "enc_distal_loose_2_tile" & t$support == 1 ] = 4
c[ t$dataset == "enc_distal_loose_2_tile" & t$support > 1 ] = 5
c[ t$dataset == "enc_ar" ] = 6
c[ t$dataset == "wang_preCRMs_tile" ] = 7
c[ t$dataset == "put_enc_proms_tile" ] = 8

t = cbind( t, c )

t = t[ t$c <= 5, ]

# Quantile filter
q = .90
filtered = t[ !( ( t$seq_gc < quantile( t$seq_gc, q ) ) & ( t$phastCons < quantile( t$phastCons, q ) ) & ( t$F < quantile( t$F, q ) ) ), ]

# Palette
## palette( c( rgb( 0.66, 0.66, 1 ), rgb( 0.33, 0.33, 1 ), rgb( 0, 0, 1 ), rgb( 1, 0.33, 0.33 ), rgb( 1, 0, 0 ), "orange", "green" ) )
## pch = c( 16, 16, 16, 16, 16, 1, 1 )

F = t$F

high_seq_gc = t[ ( t$seq_gc > quantile( t$seq_gc, q ) ) & ( t$phastCons < quantile( t$phastCons, q ) ) & ( F < quantile( F, q ) ), ]
high_phastCons = t[ ( t$seq_gc < quantile( t$seq_gc, q ) ) & ( t$phastCons > quantile( t$phastCons, q ) ) & ( F < quantile( F, q ) ), ]
high_F = t[ ( t$seq_gc < quantile( t$seq_gc, q ) ) & ( t$phastCons < quantile( t$phastCons, q ) ) & ( F > quantile( F, q ) ), ]
high_seq_gc_phastCons = t[ ( t$seq_gc > quantile( t$seq_gc, q ) ) & ( t$phastCons > quantile( t$phastCons, q ) ) & ( F < quantile( F, q ) ), ]
high_seq_gc_F = t[ ( t$seq_gc > quantile( t$seq_gc, q ) ) & ( t$phastCons > quantile( t$phastCons, q ) ) & ( F > quantile( F, q ) ), ]
high_phastCons_F = t[ ( t$seq_gc < quantile( t$seq_gc, q ) ) & ( t$phastCons > quantile( t$phastCons, q ) ) & ( F > quantile( F, q ) ), ]
high_all = t[ ( t$seq_gc > quantile( t$seq_gc, q ) ) & ( t$phastCons > quantile( t$phastCons, q ) ) & ( F > quantile( F, q ) ), ]

write.table( high_seq_gc, "F2.high_seq_gc.dat", row.names=FALSE, quote=FALSE, sep="\t" )
write.table( high_phastCons, "F2.high_phastCons.dat", row.names=FALSE, quote=FALSE, sep="\t" )
write.table( high_F, "F2.high_F.dat", row.names=FALSE, quote=FALSE, sep="\t" )
