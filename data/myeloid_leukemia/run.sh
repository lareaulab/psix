DATA='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/Giustacchini'

python ~/psix/utils/psix.py -psi $DATA/skipped_exons_psi.tab -mrna $DATA/mrna_per_event.tab -rd $DATA/rd_umap5.tab -o psix_myeloid_leukemia -k 100
