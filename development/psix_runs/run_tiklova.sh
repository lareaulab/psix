DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

python ../psix_developer.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis -k 100 

python ../psix_developer.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis_reduced_var -k 100 -pv 0.5