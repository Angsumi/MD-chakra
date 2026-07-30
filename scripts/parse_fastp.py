import json
import glob

print('Sample\tClean data (bp)\tPair-end reads\tQ30 (%)\tGC content (%)')
for f in sorted(glob.glob('fastp_qc/*_fastp_report.json')):
    with open(f) as jf:
        data = json.load(jf)
        sample = f.split('/')[-1].split('_')[0]
        summary = data['summary']['after_filtering']
        clean_bp = summary['total_bases']
        read_pairs = summary['total_reads'] // 2
        q30 = summary['q30_rate'] * 100
        gc = summary['gc_content'] * 100
        print(f'{sample}\t{clean_bp}\t{read_pairs}\t{q30:.2f}%\t{gc:.2f}%')
