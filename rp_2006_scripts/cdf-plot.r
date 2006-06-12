library(stepfun)

model 		= commandArgs()[2]
out_fname   = commandArgs()[3]
#title 		= commandArgs()[4]

pdf( out_fname, width=4, height=4, pointsize=10 )
par( mar=c(4,4,1,1) )

datasets = list( c( "Reg. training set",        "reg" ),
                 #c( "Encode promoters training set", "pos_enc_proms_chopped" ),
           		 #c( "Encode promoters (positive)", 	  "pos_enc_proms_tile" ),
           		 #c( "Encode promoters (ubiquitous)",  "ubiq_enc_proms_tile" ),
           		 #c( "Encode promoters (negative)",    "neg_enc_proms_tile" ),
           		 #c( "Encode distal training set",     "enc_distal_train_chopped" ),
           		 #c( "Encode distal",                  "enc_distal_tile" ),
           		 #c( "Encode distal (loose)",          "enc_distal_loose_tile" ),
		   		 #c( "Wang et al. validated preCRMs",  "wang_preCRMs_tile" ),
           		 c( "Exons", 						  "exons" ),
           		 c( "Bulk", 						  "all" ),
           		 c( "AR training set",                "ar" ) )
           		 #c( "ENCODE AR training set",         "enc_ar" ),
		   		 #c( "REG -- tss",                     "reg_tss_tile" ),
 				 #c( "REG -- proximal",                "reg_proximal_tile" ),
                 #c( "REG -- genic",                   "reg_genic_tile" ),
				 #c( "REG -- distal",                  "reg_distal_tile" ) )

data = list()
descs = NULL

for ( l in datasets )
{
    desc = l[1]; name = l[2]
    d = scan( paste( name, ".", model, ".scores", sep="" ) )
    d = d[ !is.na( d ) ]
    data[[name]] = d
    descs = c( descs, desc )
}

ndatasets = length( datasets )           

lims = range( data["reg"], data["ar"] ) 
#lims = c( -.1, .2 )
#dummy = seq( lims[1], lims[2], 0.001 )
print( lims )
dummy = seq( lims[1], lims[2], (lims[2]-lims[1])/1000 )

#cols=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
cols = rainbow( ndatasets + 1, 1.0, 0.7 )
ltys=rep( 1, ndatasets )

#cols=c("black","gray","gray","black")
#ltys=c(1,1,2,2)

plot( dummy, xlim=lims, ylim=range(0,1), xlab="Score", ylab="Cumulative Distribution" )
# lines( dummy, ecdf(reg_shuff)(dummy), lwd=2, col="magenta" )
# lines( dummy, ecdf(enhancers)(dummy), lwd=1, col="darkgoldenrod" )

color = 0

for ( l in datasets )
{
    desc = l[1]; name = l[2]
    d = data[[name]]
    lines( dummy, ecdf(d)(dummy), lwd=2, lty=1, col=cols[ ( color = color + 1 ) ] )
}

legend( "bottomright", NULL, descs, col=cols, lwd=2, lty=1, bty="n", x.intersp=1.5, y.intersp=1.5 )

n = NULL
v = NULL

for ( i in (1:ndatasets) )
{
	v = c( data[[i]], v )
	n = c( rep( datasets[[i]][[1]], length( data[[i]] ) ), n )
}

par( mar=c(4,16,1,1), las=1 )
boxplot( v ~ n, list( n=as.factor(n), v=v ), horizontal=T, varwidth=T )

dev.off()
