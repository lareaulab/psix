DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

python ../utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab -o tiklova_neurogenesis -k 100 -a 10

#python ../utils/psix_cross.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc2.tab  -p tiklova_neurogenesis -o tiklova_neurogenesis -k 100 -a 10



python ../utils/psix.py -psi $DATA/chen/skipped_exons_psi.tab -mrna $DATA/chen/mrna_per_event.tab -rd $DATA/chen/rd_pc2.tab -o chen -k 30 -a 10

#python ../utils/psix_cross.py -psi $DATA/chen/skipped_exons_psi.tab -mrna $DATA/chen/mrna_per_event.tab -rd $DATA/chen/rd_pc2.tab  -p chen -o chen -k 30 -a 10

####

python ../utils/psix.py -psi $DATA/tiklova/skipped_exons_psi.tab -mrna $DATA/tiklova/mrna_per_event.tab -rd $DATA/tiklova/rd_pc2.tab -o tiklova -k 100 -a 10

#python ../utils/psix_cross.py -psi $DATA/tiklova/skipped_exons_psi.tab -mrna $DATA/tiklova/mrna_per_event.tab -rd $DATA/tiklova/rd_pc2.tab  -p tiklova -o tiklova -k 100 -a 10





python ../utils/psix.py -psi $DATA/song/skipped_exons_psi.tab -mrna $DATA/song/mrna_per_event.tab -rd $DATA/song/rd_pc2.tab -o song -k 30 -a 10

#python ../utils/psix_cross.py -psi $DATA/song/skipped_exons_psi.tab -mrna $DATA/song/mrna_per_event.tab -rd $DATA/song/rd_pc2.tab  -p song -o song -k 30 -a 10




python ../utils/psix.py -psi $DATA/tiklova_neurogenesis/skipped_exons_psi.tab -mrna $DATA/tiklova_neurogenesis/mrna_per_event.tab -rd $DATA/tiklova_neurogenesis/rd_pc5.tab -o tiklova_neurogenesis_pc5 -k 30 -a 10


python ../utils/psix.py -psi $DATA/tiklova/skipped_exons_psi.tab -mrna $DATA/tiklova/mrna_per_event.tab -rd $DATA/tiklova/rd_pc5.tab -o tiklova_pc5 -k 30 -a 10



