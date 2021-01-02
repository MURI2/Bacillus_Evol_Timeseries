#############################
#
# This script processes output from cluster to obtain final set
# of mutation trajectories that are used for downstream analysis
#
#############################

import os, imp, sys
from datetime import date
import phylo_tools as pt


likelihood_directory = pt.get_path() + '/data/timecourse_likelihood/'
filter_snps_path = pt.get_path() + '/Python/filter_snps_and_calculate_depth.py'
pvals_python_path = pt.get_path() + '/Python/annotate_pvalues.py'
filter_pvals_python_path = pt.get_path() + '/Python/filter_pvalues.py'

combined_annotate_timecourse_path = pt.get_path() + '/Python/combined_annotate_timecourse.py'

well_mixed_hmm_wrapper = pt.get_path() + '/Python/calculate_well_mixed_hmm_wrapper.py'


treatments = ['0', '1', '2']
reps = ['1', '2', '3', '4', '5']
#treatments = ['0']
#strains = ['B', 'C', 'D', 'F', 'J', 'P']
# fix S strain
strains = ['S']

for strain in strains:
    sys.stdout.write("\nProcessing Strain %s...\n" % (strain))
    for treatment in treatments:
        for rep in reps:

            print('%s%s%s' % (treatment, strain, rep))

            if treatment+strain+rep in pt.populations_to_ignore:
                continue

            merged_timecourse_filename = '%s%s%s_merged_timecourse.bz' % (treatment, strain, rep)
            depth_timecourse_filename = '%s%s%s_depth_timecourse.bz' % (treatment, strain, rep)
            snp_timecourse_filename = '%s%s%s_snp_timecourse.bz' % (treatment, strain, rep)
            indel_timecourse_filename = '%s%s%s_indel_timecourse.bz' % (treatment, strain, rep)
            likelihood_timecourse_filename = '%s%s%s_likelihood_timecourse.bz' % (treatment, strain, rep)
            #likelihood_timecourse_clean_filename = '%s%s%s_likelihood_timecourse_clean.bz' % (treatment, strain, rep)
            annotated_timecourse_filename = '%s%s%s_annotated_timecourse.bz' % (treatment, strain, rep)

            merged_timecourse_path = pt.get_path() + '/data/timecourse_merged/' + merged_timecourse_filename
            depth_timecourse_path = pt.get_path() + '/data/timecourse_depth/' + depth_timecourse_filename
            snp_timecourse_path = pt.get_path() + '/data/timecourse_snp/' + snp_timecourse_filename
            indel_timecourse_path = pt.get_path() + '/data/timecourse_indel/' + indel_timecourse_filename
            likelihood_timecourse_path = pt.get_path() + '/data/timecourse_likelihood/' + likelihood_timecourse_filename
            #likelihood_timecourse_clean_path = pt.get_path() + '/data/timecourse_likelihood_clean/' + likelihood_timecourse_clean_filename
            annotated_timecourse_path = pt.get_path() + '/data/timecourse_final/' + annotated_timecourse_filename

            # Filter SNPs and calculate avg depth per sample
            #sys.stdout.write('Filtering SNPS and calculating depth...\n')
            #return_val = os.system('python %s %s %s %s %s' % (filter_snps_path, merged_timecourse_path, depth_timecourse_path, snp_timecourse_path, strain))
            #if return_val==0:
            #    sys.stdout.write('Done!\n')
            #else:
            #    sys.stdout.write("Error!\n")


            # Call indels
            #sys.stdout.write("Calling indels...\n")
            #call_indels_path = pt.get_path() + '/Python/call_indels.py'
            #return_val = os.system('python %s %s %s' % (call_indels_path, merged_timecourse_path, indel_timecourse_path))
            #if return_val==0:
            #    sys.stdout.write('Done!\n')
            #else:
            #    sys.stdout.write("Error!\n")


            # Get pvalues
            #sys.stdout.write("Calculating pvalues...\n")
            #return_val = os.system('python %s %s %s %s %s' % \
            #    (pvals_python_path, depth_timecourse_path, snp_timecourse_path, indel_timecourse_path, likelihood_timecourse_path))
            #if return_val==0:
            #    sys.stdout.write('Done!\n')
            #else:
            #    sys.stdout.write("Error!\n")



    sys.stdout.write("\n\nTrajectory post-processing output for energy-limited evolution sequencing project\n")
    sys.stdout.write("Date: %s\n\n" % str(date.today()))

    # Filter and annotate timecourse
    sys.stdout.write('Performing final filtering and annotation step...\n')
    return_val = os.system('python %s %s' % (combined_annotate_timecourse_path, strain))
    if return_val==0:
        sys.stdout.write('Done!\n')
    else:
        sys.stdout.write("Error!\n")

    # Infer trajectory states in well-mixed HMM
    sys.stdout.write('Inferring trajectory states in well-mixed HMM...\n')
    return_val = os.system('python %s %s' % (well_mixed_hmm_wrapper, strain))
    if return_val==0:
        sys.stdout.write('Done!\n')
    else:
        sys.stdout.write("Error!\n")




# Position, Gene, Allele, Annotation, Test statistic, P-value,
# Deletion index, Fold reduction, Deletion P-value, Duplication index,
# Fold increase, Duplication pvalue, Passed?, AC:0, D
