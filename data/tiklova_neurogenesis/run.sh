DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

#python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis -k 100 -a 10

python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis_extended -k 100 -a 10 -n 0.1

python ~/psix/utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc5.tab -o tiklova_neurogenesis_pc5 -k 100 -a 10


#python ../utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis_extended -k 100 -a 10 -n 0.1
