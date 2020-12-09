DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

python ~/psix/utils/psix.py -psi tiklova_random_subsampled_psi.tab -mrna tiklova_random_subsampled_mrna.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_random_subsampled -k 100 -a 10 -n 0.01
 
 python ~/psix/utils/psix.py -psi tiklova_random_psi.tab -mrna tiklova_random_mrna.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_random -k 100 -a 10 -n 0.01
