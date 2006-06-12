scores_fname = commandArgs()[2]
out_fname    = commandArgs()[3]

source( "~dcking/work/calibration/binary_classification/Rutil/classifiers.R" )

lineROC = function(dl,...) {
    lines(x= dl$FP / ( dl$TN + dl$FP),       # specificity
         y= dl$TP / ( dl$TP + dl$FN),       # sensitivity
        lwd=2, ... )
    perf = performance(dl)
    # points( 1-perf$Sp, perf$Sn, pch=8, ... )
}

nbreaks = 2048

data_mcs = read.table( "~dcking/work/calibration/binary_classification/evaluation/mcs.hmr.lab" )
data_mcs = data_mcs[ data_mcs[,2] > 0, ] # Filter so that only scores > 0 are considered
table_mcs = discriminate( data_mcs, nbreaks )

data_phast = read.table( "~dcking/work/calibration/binary_classification/evaluation/phc.hmr.lab" )
table_phast = discriminate( data_phast, nbreaks )

data_rp_hmr = read.table( "~dcking/work/calibration/binary_classification/evaluation/rp.hmr.lab" )
table_rp_hmr = discriminate( data_rp_hmr, nbreaks )

# data_rp_hmd = read.table( "~dcking/work/calibration/binary_classification/evaluation/rp.hg17.mm5.canFam1.lab" )
# table_rp_hmd = discriminate( data_rp_hmd, nbreaks )

#data_rp_5way = read.table( "~dcking/work/calibration/binary_classification/evaluation/5species_hbb.full_scores.lab" )
#data_rp_5way = read.table( "~/work/rp_training_data/hg17/5way_new/hbbc.hc_125_1_best.window_scores.lab" )
#data_rp_5way = read.table( "~/work/rp_training_data/hg17/5way_new/hbbc.knn_80_3_best.window_scores.lab" )
#data_rp_5way = read.table( "~/work/rp_training_data/hg17/5way_new/hbbc.knn_80_4_best.window_scores.lab" )
data_rp_5way = read.table( "~/work/rp_training_data/hg17/5way_new/hbbc.knn_ent_agglom_best.window_scores.lab" )
#data_rp_5way = read.table( "~/work/rp_training_data/hg17/5way_new/hbbc.knn_80_inti_best.window_scores.lab" )

data_rp_5way = read.table( scores_fname )

table_rp_5way = discriminate( data_rp_5way, nbreaks )

# Plotting

pdf( out_fname, width=4, height=4, pointsize=10 )
par( mar=c(4,4,1,1) )

plot( 0, 0, xlim=c(0,1), ylim=c(0,1), xlab="1 - Specificity", ylab="Sensitivity", type="n" )
abline( 0, 1, lty=3 )

cols=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
ltys=c(1,1,1,1)

#cols=c("black","gray","gray","black")
#ltys=c(2,1,2,1)

lineROC( table_mcs, col=cols[1], lty=ltys[1] )
lineROC( table_phast, col=cols[2], lty=ltys[2])
lineROC( table_rp_hmr, col=cols[3], lty=ltys[3])
# lineROC( table_rp_hmd, col=5 )
lineROC( table_rp_5way, col=cols[4], lty=ltys[4])

legend( "bottomright", NULL, c( "MCS", "PhastCons", "Human/Mouse/Rat RP", "Five Species RP" ), col=cols, lwd=2, lty=ltys, bty="n", x.intersp=1.5, y.intersp=1.5 )

