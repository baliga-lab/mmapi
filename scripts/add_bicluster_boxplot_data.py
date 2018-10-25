#!/usr/bin/env python3
import sys
import MySQLdb
import numpy as np
import pandas

"""
For box plot data, we need to extract min, max, lower quartile, upper quartile and median
from the gene expression matrix
"""

if __name__ == '__main__':
    if len(sys.argv) > 1:
        db = MySQLdb.connect(host="localhost",
                             user="egrin2",
                             passwd="egrin2",
                             db="mm_api")
        cursor = db.cursor()
        df = pandas.read_csv(sys.argv[1], sep='\t', index_col=0, header=0)
        # now for each patient, get all the values of the genes in the bicluster
        # meaning, for each bicluster, extract the gene list, and then iterate
        # over the patients to get all the expression values in the cluster
        try:
            cursor.execute('select id,name from patients')
            all_patients = [(row[0], row[1]) for row in cursor.fetchall()]
            cursor.execute('select id,name from biclusters')
            for bc_id, bc_name in cursor.fetchall():
                cursor2 = db.cursor()
                cursor2.execute('select g.ensembl_id from bicluster_genes bg join genes g on bg.gene_id=g.id where bg.bicluster_id=%s', [bc_id])
                bc_genes = [row[0] for row in cursor2.fetchall()]

                for patient_id, patient_name in all_patients:
                    values = df.loc[bc_genes][patient_name]
                    median_value = np.median(values)
                    min_value = np.min(values)
                    max_value = np.max(values)
                    lower_quartile = np.percentile(values, 25)
                    upper_quartile = np.percentile(values, 75)
                    cursor2.execute('insert into bicluster_boxplot_data (bicluster_id,patient_id,median,min_value,max_value,lower_quartile,upper_quartile) values (%s,%s,%s,%s,%s,%s,%s)',
                                    [bc_id, patient_id, median_value, min_value, max_value,
                                     lower_quartile, upper_quartile])
            db.commit()

            print("done")
        finally:
            cursor.close()
            db.close()
    else:
        print('please provide gene expression file')
