
import os
import re
import subprocess

def copy_gm_key():
    if not os.path.exists('./.gm_key'):
        homedirkey = os.path.expanduser('~') + '/.gm_key'
        if not os.path.exists(homedirkey):
            print('GeneMark license not found in home directory. Please download .gm_key from https://genemark.bme.gatech.edu/license_download.cgi')
        else:
            subprocess.call(['cp',homedirkey,'./.gm_key'])
            print('Copied GeneMark license to working directory')
    else:
        print('GeneMark license found in funQCD directory')

def make_samples_csv(short_read_metadata,long_read_metadata,short_read_pass_dir,long_read_pass_dir):
    pairs = find_pair(short_read_metadata,long_read_metadata,short_read_pass_dir,long_read_pass_dir)
    # this should contain pairs of short read and long read names in this format:
    # short_internal,long_internal,short_isolate,long_isolate
    if os.path.exists('config/samples.csv'):
        print('Overwriting config/samples.csv with new version')
    with open('config/samples.csv','w') as fhout:
        _ = fhout.write('ShortRead_sample_id,LongRead_sample_id,ShortRead_isolate_name,LongRead_isolate_name\n')
        for pairdata in pairs:
            _ = fhout.write(','.join(pairdata) + '\n')

def find_pair(short_read_metadata,long_read_metadata,short_read_pass_dir,long_read_pass_dir):
    pairs = []
    short_read_dict = make_dict(short_read_metadata)
    long_read_dict = make_dict(long_read_metadata)
    pass_ShortRead_InternalName = set([re.split('_R[1,2].fastq.gz',x)[0] for x in os.listdir(short_read_pass_dir)])
    # this is a set of each name present in the directory of passing short read files
    # these are the internal names, such as UM_Caur_10
    pass_LongRead_InternalName = [x.split('.fastq.gz')[0] for x in os.listdir(long_read_pass_dir)]
    # this is a list of the internal names of the samples that passed long read assembly
    pass_LongRead_ISOLATEname = [long_read_dict.get(x) for x in pass_LongRead_InternalName]
    if None in pass_LongRead_ISOLATEname:
        print(f'The files present in {long_read_pass_dir} were not found in {long_read_metadata}!')
        quit(1)
    # this converts each long read internal name into the corresponding isolate name
    # both lists should have identical indexing (InternalName[5] corresponds to ISOLATEname[5])
    for short_IntName in pass_ShortRead_InternalName:
        data_pair = None
        if short_IntName not in short_read_dict:
            print(f'The files present in {short_read_pass_dir} were not found in {short_read_metadata}!')
            quit(1)
        short_ISOLATE_NAME = short_read_dict[short_IntName]
        # this converts the short read internal name to the isolate name
        for ind in range(len(pass_LongRead_ISOLATEname)):
            # for each index in the list of long read isolate names:
            long_IntName = pass_LongRead_InternalName[ind]
            long_ISOLATE_NAME = pass_LongRead_ISOLATEname[ind]
            if long_ISOLATE_NAME == short_ISOLATE_NAME:
                # if the long read isolate name for that index is the same as the short read isolate name:
                data_pair = (short_IntName,long_IntName,short_ISOLATE_NAME,long_ISOLATE_NAME)
                # then pair the short read internal name with the long read internal name of the same index
                # (it's OK if there are multiple matches, as data_pair will just be overwritten)
        if data_pair is not None:
            pairs.append(data_pair)
    # this should contain pairs of internal names whose isolate names match
    return(pairs)


def make_dict(infile):
    d = {}
    with open(infile,'r') as fh:
        next(fh)
        for line in fh:
            # the first two entries in each line should be the sample ID (such as CaTO_703), followed by the internal name (UM_Caur_1)
            # the same sample ID does not always have the same internal name assigned to it
            # for example, we could have CaTO_713=UM_Caur_10 for Illumina data and CaTO_713=UM_Caur_nano_2 for Nanopore data
            # internal names should all be unique, while sample IDs are not
            sample_id,internal_name = line.strip().split('\t')[:2]
            d[internal_name] = sample_id
    return(d)

if __name__ == "__main__":
    short_read_metadata = '/nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/metadata/sample_lookup/UM_illumina_sample_lookup.txt'
    long_read_metadata = '/nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/metadata/sample_lookup/UM_ONT_sample_lookup.txt'
    short_read_batch_name = '2025_01_29_UM_illumina'
    long_read_batch_name = '2025_01_29_UM_ONT'
    # only the text above needs editing
    short_read_pass_dir = '/nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/illumina_fastq/' + short_read_batch_name + '/passed_qc_samples/'
    long_read_pass_dir = '/nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/ONT/' + long_read_batch_name + '/passed_qc_samples/'
    copy_gm_key()
    make_samples_csv(short_read_metadata,long_read_metadata,short_read_pass_dir,long_read_pass_dir)
    print(f'Remember to add {short_read_pass_dir} and {long_read_pass_dir} to config/config.yaml!')
