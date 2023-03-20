import pathlib
import pandas
import os

SAMPLES = None
MAP_SAMPLE_TO_INPUT = None


def process_sample_sheet():
    # sample_sheet_file = pathlib.Path('config/samples.tsv').resolve(strict=True)
    sample_sheet_file = pathlib.Path(config["samples"]).resolve(strict=True)
    sample_sheet = pandas.read_csv(
        sample_sheet_file,
        sep="\t",
        header=0
    )

    if sample_sheet['sample'].duplicated().any():
        raise ValueError('Duplicated entries in "sample" column')

    for col in ['gfa', 'strandseq_dir']:
        sample_sheet[col] = sample_sheet[col].map(pathlib.Path)

    samples = set()
    sample_input = dict()

    for row in sample_sheet.itertuples(index=False):
        ss_dir = row.strandseq_dir.resolve(strict=True)
        ss_suffix = row.strandseq_suffix
        ss_files = set(ss_dir.glob(f"**/*{ss_suffix}"))

        ss_files = [str(f) for f in sorted(ss_files) if f.is_file()]
        if len(ss_files) < 1:
            raise ValueError(f"No files found underneath {ss_dir}")

        ss_libs = extract_libs(ss_dir, ss_suffix)

        if len(ss_libs) *  2 != len(ss_files):
            print(row.sample)
            print(sorted(ss_libs))
            print(len(ss_libs))

        assert(len(ss_libs) *  2 == len(ss_files)), 'Incongruent number of Strand-seq libs and files'
        samples.add(row.sample)
        sample_input[row.sample] = dict(
            gfa=str(row.gfa.resolve(strict=True)),
            strandseq_dir=str(ss_dir),
            strandseq_suffix=ss_suffix,
            strandseq=ss_files,
            strandseq_libs=ss_libs
        )

    global SAMPLES
    SAMPLES = sorted(samples)
    global MAP_SAMPLE_TO_INPUT
    MAP_SAMPLE_TO_INPUT = sample_input

def extract_libs(ss_dir, ss_suffix):
    libs, = glob_wildcards(pathlib.Path(ss_dir,"{lib}_1" + ss_suffix))
    if not libs:
        raise FileNotFoundError(f'No Strand-Seq files found in {ss_dir}. Strand-Seq files are expected to have filenames of the form: {pathlib.Path(ss_dir,"{lib}_{1,2}" + ss_suffix)}')
    return libs


process_sample_sheet()