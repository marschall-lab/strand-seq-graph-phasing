import pathlib
import pandas

SAMPLES = None
MAP_SAMPLE_TO_INPUT = None


def process_sample_sheet():
    # sample_sheet_file = pathlib.Path('config/samples-hgsvc-verkko14.tsv').resolve(strict=True)
    sample_sheet_file = pathlib.Path(config["samples"]).resolve(strict=True)
    sample_sheet = pandas.read_csv(
        sample_sheet_file,
        sep="\t",
        header=0
    )

    # TODO add schema checks for file

    if sample_sheet['sample'].duplicated().any():
        raise ValueError('Duplicated entries in "sample" column')

    samples = set()
    sample_input = dict()

    for row in sample_sheet.itertuples(index=False):
        print('Processing sample: ' + row.sample)

        samples.add(row.sample)
        sample_input[row.sample] = dict()

       # Strand-seq files
        # ss_dir = pathlib.Path('/Users/henglinm/Documents/wd/GM19317').resolve(strict=True)
        ss_dir = pathlib.Path(row.strandseq_dir).resolve(strict=True)
        ss_files = set(ss_dir.glob(f"**/*.gz"))

        ss_files = [str(f) for f in sorted(ss_files) if f.is_file()]
        if len(ss_files) < 1:
            raise ValueError(f"No files found underneath {ss_dir}")

        sseq_pairs = organize_sseq_files(ss_files)

        sample_input[row.sample]['strandseq_libs'] = list(sseq_pairs.keys())
        sample_input[row.sample]['strandseq_pairs'] = sseq_pairs

        # GFA
        sample_input[row.sample]['gfa'] = str(pathlib.Path(row.gfa).resolve(strict=True))

        # hpc
        sample_input[row.sample]['hpc'] = row.hpc

        # assembler
        sample_input[row.sample]['assembler'] = row.assembler

        # assembler
        sample_input[row.sample]['cluster_PAR_with_haploid'] = row.cluster_PAR_with_haploid


        # Coverage ~ Optional
        if not pandas.isna(row.coverage):
            sample_input[row.sample]['coverage'] = str(pathlib.Path(row.coverage).resolve(strict=True))
        elif row.assembler == 'hifiasm':
            sample_input[row.sample]['coverage'] = f'hifiasm_hifi_coverage/{row.sample}_hifiasm_hifi_coverage.tsv'


    global SAMPLES
    SAMPLES = sorted(samples)
    global MAP_SAMPLE_TO_INPUT
    MAP_SAMPLE_TO_INPUT = sample_input


def is_different_in_one_location(s1, s2):
    """
    Returns True if two strings are different only in one location
    and the characters at the differing location are 1 or 2.

    Parameters:
    s1 (str): The first string.
    s2 (str): The second string.

    Returns:
    bool: True if the two strings are different only in one location
          and the characters at the differing location are 1 or 2,
          False otherwise.
    """
    if len(s1) != len(s2):
        return False

    allowed_values = ['1', '2']
    different_locations = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            different_locations += 1
            if different_locations > 1 or (s1[i] not in allowed_values) or (s2[i] not in allowed_values):
                return False

    return different_locations == 1

def organize_sseq_files(ss_files, id_size = 6):

    sseq_pairs = [(ss_files[i], ss_files[i + 1]) for i in range(0,len(ss_files) - 1,2)]
    if not all(is_different_in_one_location(x[0], x[1]) for x in sseq_pairs):
        raise ValueError('Strandseq files are mismatched')

    # with artificial library ids, because I am lazy right now
    sseq_lib_dict = {}
    for i in range(len(sseq_pairs)):
        id = str(i).zfill(id_size)
        sseq_lib_dict[id] = {}
        sseq_lib_dict[id]['file1'] = sseq_pairs[i][0]
        sseq_lib_dict[id]['file2'] = sseq_pairs[i][1]

    return sseq_lib_dict

process_sample_sheet()