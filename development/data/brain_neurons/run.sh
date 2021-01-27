DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation'

#python ~/psix/utils/psix.py -psi $DATA/brain_neurons/skipped_exons_psi.tab -mrna $DATA/brain_neurons/mrna_per_event.tab -rd $DATA/brain_neurons/rd_pc30.tab -o brain_neurons -k 100 -a 10

python ~/psix/utils/psix.py -psi $DATA/brain_neurons/skipped_exons_psi.tab -mrna $DATA/brain_neurons/mrna_per_event.tab -rd $DATA/brain_neurons/tm_scvi_rd.tab -o brain_neurons_scvi -k 100
