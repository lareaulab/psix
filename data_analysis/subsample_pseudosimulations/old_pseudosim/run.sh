DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

python ~/psix/utils/psix.py -psi tiklova_subsampled_exon_psi.tab.gz -mrna tiklova_subsampled_exon_mrna.tab.gz -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis -k 100 -a 10
