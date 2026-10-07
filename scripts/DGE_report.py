#!/usr/bin/env python3
"""Summarise collated DE/GSEA comparison outputs into interpretation tables and HTML reports."""
import argparse, html, json, re
from pathlib import Path
import numpy as np
import pandas as pd


def safe_filename(x):
    x = re.sub(r"[^\w.-]+", "_", str(x).strip())
    return re.sub(r"_+", "_", x).strip("_.")


def args_parser():
    p=argparse.ArgumentParser()
    p.add_argument('--comparison_dir', required=True, help='Output directory made by compare_DGE.py')
    p.add_argument('--reference_analysis', required=True)
    p.add_argument('--analyses', nargs='+', required=True)
    p.add_argument('--report_name', default='interpretation_report')
    p.add_argument('--top_n', type=int, default=25)
    p.add_argument('--min_abs_log2fc', type=float, default=0.5)
    p.add_argument('--min_abs_nes', type=float, default=1.0)
    p.add_argument('--fdr', type=float, default=0.05)
    return p.parse_args()


def read_csv(path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def numeric(s): return pd.to_numeric(s, errors='coerce')


def median_abs_dev_from_ref(df, analyses, ref, metric):
    refcol=f'{ref}_{metric}'
    diffs=[]
    for a in analyses:
        c=f'{a}_{metric}'
        if a!=ref and c in df.columns and refcol in df.columns:
            diffs.append((numeric(df[c])-numeric(df[refcol])).abs())
    return pd.concat(diffs, axis=1).median(axis=1, skipna=True) if diffs else pd.Series(np.nan,index=df.index)


def classify_gene_table(df, analyses, ref, fdr, min_fc):
    out=df.copy()
    lfc=[f'{a}_log2FoldChange' for a in analyses if f'{a}_log2FoldChange' in out]
    sig=[f'{a}_FDR' for a in analyses if f'{a}_FDR' in out]
    out['n_available_effects']=out[lfc].notna().sum(axis=1)
    out['n_significant']=sum(numeric(out[c]).lt(fdr) for c in sig) if sig else 0
    if lfc:
        signs=np.sign(out[lfc].apply(numeric))
        out['n_positive']=signs.gt(0).sum(axis=1); out['n_negative']=signs.lt(0).sum(axis=1)
        out['mean_abs_log2FC']=out[lfc].apply(numeric).abs().mean(axis=1)
        out['effect_range']=out[lfc].apply(numeric).max(axis=1)-out[lfc].apply(numeric).min(axis=1)
    out['median_abs_delta_from_reference']=median_abs_dev_from_ref(out,analyses,ref,'log2FoldChange')
    rlfc=f'{ref}_log2FoldChange'; rfdr=f'{ref}_FDR'
    out['reference_significant']=numeric(out[rfdr]).lt(fdr) if rfdr in out else False
    out['reference_large_effect']=numeric(out[rlfc]).abs().ge(min_fc) if rlfc in out else False
    out['agreement_priority']=(
        out['reference_significant'].astype(int)*100 +
        out.get('all_same_direction',False).fillna(False).astype(int)*25 +
        out['n_available_effects']*5 + out['n_significant']*5 -
        out['median_abs_delta_from_reference'].fillna(10)
    )
    out['standout_type']=np.select([
        out['reference_significant'] & out.get('all_same_direction',False).fillna(False),
        out['reference_significant'] & ~out.get('all_same_direction',False).fillna(False),
        (~out['reference_significant']) & (out['n_significant']>0),
        out.get('direction_class',pd.Series('',index=out.index)).eq('MIXED_DIRECTION')
    ],['REFERENCE_HIT_CONCORDANT','REFERENCE_HIT_DISCORDANT','SUBSET_SPECIFIC_HIT','DIRECTION_DISCORDANT'],default='OTHER')
    return out


def classify_pathway_table(df, analyses, ref, fdr, min_nes):
    out=df.copy(); nes=[f'{a}_NES' for a in analyses if f'{a}_NES' in out]; sig=[f'{a}_FDR q-val' for a in analyses if f'{a}_FDR q-val' in out]
    out['n_available_NES']=out[nes].notna().sum(axis=1)
    out['n_significant']=sum(numeric(out[c]).lt(fdr) for c in sig) if sig else 0
    if nes:
        N=out[nes].apply(numeric); signs=np.sign(N)
        out['n_positive']=signs.gt(0).sum(axis=1); out['n_negative']=signs.lt(0).sum(axis=1)
        out['mean_abs_NES']=N.abs().mean(axis=1); out['NES_range']=N.max(axis=1)-N.min(axis=1)
    out['median_abs_delta_from_reference']=median_abs_dev_from_ref(out,analyses,ref,'NES')
    rn=f'{ref}_NES'; rf=f'{ref}_FDR q-val'
    out['reference_significant']=numeric(out[rf]).lt(fdr) if rf in out else False
    out['reference_strong']=numeric(out[rn]).abs().ge(min_nes) if rn in out else False
    out['agreement_priority']=(out['reference_significant'].astype(int)*100 + out.get('all_same_direction',False).fillna(False).astype(int)*25 + out['n_available_NES']*5 + out['n_significant']*5 - out['median_abs_delta_from_reference'].fillna(10))
    out['standout_type']=np.select([
        out['reference_significant'] & out.get('all_same_direction',False).fillna(False),
        out['reference_significant'] & ~out.get('all_same_direction',False).fillna(False),
        (~out['reference_significant']) & (out['n_significant']>0),
        out.get('direction_class',pd.Series('',index=out.index)).eq('MIXED_DIRECTION')
    ],['REFERENCE_HIT_CONCORDANT','REFERENCE_HIT_DISCORDANT','SUBSET_SPECIFIC_HIT','DIRECTION_DISCORDANT'],default='OTHER')
    return out


def html_table(df, cols, n=15):
    cols=[c for c in cols if c in df.columns]
    x=df[cols].head(n).copy()
    for c in x.columns:
        if pd.api.types.is_numeric_dtype(x[c]): x[c]=x[c].map(lambda v:'' if pd.isna(v) else f'{v:.3g}')
    return x.to_html(index=False, escape=True, border=0, classes='data')


def main():
    a=args_parser(); root=Path(a.comparison_dir); out=root/'reports'/safe_filename(a.report_name); out.mkdir(parents=True,exist_ok=True)
    manifest={'comparison_dir':str(root.resolve()),'reference_analysis':a.reference_analysis,'analyses':a.analyses,'top_n':a.top_n,'fdr':a.fdr,'min_abs_log2fc':a.min_abs_log2fc,'min_abs_nes':a.min_abs_nes}
    (out/'report_manifest.json').write_text(json.dumps(manifest,indent=2))

    de_dir=root/'DE'/'by_cell_type'; gsea_dir=root/'GSEA'/'by_cell_type'
    gene_rows=[]; pathway_rows=[]; sections=[]

    for f in sorted(de_dir.glob('*.csv')):
        cell=f.stem; df=classify_gene_table(pd.read_csv(f),a.analyses,a.reference_analysis,a.fdr,a.min_abs_log2fc)
        df.insert(0,'cell_type',cell); gene_rows.append(df)
        conc=df[df['standout_type']=='REFERENCE_HIT_CONCORDANT'].sort_values(['agreement_priority'],ascending=False)
        disc=df[df['standout_type'].isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT'])].sort_values(['agreement_priority'],ascending=False)
        subset=df[df['standout_type']=='SUBSET_SPECIFIC_HIT'].sort_values(['n_significant','mean_abs_log2FC'],ascending=False)
        sections.append((cell,conc,disc,subset))
    genes=pd.concat(gene_rows,ignore_index=True,sort=False) if gene_rows else pd.DataFrame()
    if not genes.empty:
        genes.to_csv(out/'gene_interpretation_master.csv',index=False)
        genes[genes['standout_type']=='REFERENCE_HIT_CONCORDANT'].sort_values('agreement_priority',ascending=False).to_csv(out/'genes_reference_hits_concordant.csv',index=False)
        genes[genes['standout_type'].isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT'])].to_csv(out/'genes_direction_discordant.csv',index=False)
        genes[genes['standout_type']=='SUBSET_SPECIFIC_HIT'].sort_values(['n_significant','mean_abs_log2FC'],ascending=False).to_csv(out/'genes_subset_specific.csv',index=False)

    for cell_dir in sorted(gsea_dir.glob('*')):
        if not cell_dir.is_dir(): continue
        for f in sorted(cell_dir.glob('*.csv')):
            df=classify_pathway_table(pd.read_csv(f),a.analyses,a.reference_analysis,a.fdr,a.min_abs_nes)
            df.insert(0,'gene_set',f.stem); df.insert(0,'cell_type',cell_dir.name); pathway_rows.append(df)
    paths=pd.concat(pathway_rows,ignore_index=True,sort=False) if pathway_rows else pd.DataFrame()
    if not paths.empty:
        paths.to_csv(out/'pathway_interpretation_master.csv',index=False)
        paths[paths['standout_type']=='REFERENCE_HIT_CONCORDANT'].sort_values('agreement_priority',ascending=False).to_csv(out/'pathways_reference_hits_concordant.csv',index=False)
        paths[paths['standout_type'].isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT'])].to_csv(out/'pathways_direction_discordant.csv',index=False)
        paths[paths['standout_type']=='SUBSET_SPECIFIC_HIT'].sort_values(['n_significant','mean_abs_NES'],ascending=False).to_csv(out/'pathways_subset_specific.csv',index=False)

    # Compact per-cell summary
    summary=[]
    cells=sorted(set(genes.get('cell_type',[]))|set(paths.get('cell_type',[])))
    for cell in cells:
        g=genes[genes.cell_type==cell] if not genes.empty else pd.DataFrame(); p=paths[paths.cell_type==cell] if not paths.empty else pd.DataFrame()
        summary.append({'cell_type':cell,'reference_sig_genes':int(g.get('reference_significant',pd.Series(dtype=bool)).sum()),'concordant_reference_genes':int((g.get('standout_type',pd.Series(dtype=str))=='REFERENCE_HIT_CONCORDANT').sum()),'discordant_genes':int(g.get('standout_type',pd.Series(dtype=str)).isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT']).sum()),'subset_specific_genes':int((g.get('standout_type',pd.Series(dtype=str))=='SUBSET_SPECIFIC_HIT').sum()),'reference_sig_pathways':int(p.get('reference_significant',pd.Series(dtype=bool)).sum()),'concordant_reference_pathways':int((p.get('standout_type',pd.Series(dtype=str))=='REFERENCE_HIT_CONCORDANT').sum()),'discordant_pathways':int(p.get('standout_type',pd.Series(dtype=str)).isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT']).sum()),'subset_specific_pathways':int((p.get('standout_type',pd.Series(dtype=str))=='SUBSET_SPECIFIC_HIT').sum())})
    summary=pd.DataFrame(summary); summary.to_csv(out/'cell_type_interpretation_summary.csv',index=False)

    # HTML report
    css='''<style>body{font-family:Arial,sans-serif;color:#263238;margin:34px;max-width:1400px}h1,h2{color:#244a73}h3{color:#334e68;border-bottom:1px solid #d9e2ec;padding-bottom:4px}.note{background:#eef4f8;padding:12px;border-left:4px solid #4472a3}.data{border-collapse:collapse;font-size:12px;width:100%;margin:8px 0 22px}.data th{background:#244a73;color:white;padding:6px}.data td{border:1px solid #d9e2ec;padding:5px}.data tr:nth-child(even){background:#f7f9fb}.small{color:#66788a;font-size:12px}</style>'''
    H=[f'<html><head><meta charset="utf-8"><title>{html.escape(a.report_name)}</title>{css}</head><body>',f'<h1>{html.escape(a.report_name)}</h1>',f'<div class="note"><b>Reference analysis:</b> {html.escape(a.reference_analysis)}<br><b>Compared analyses:</b> {html.escape(", ".join(a.analyses))}<br>Agreement is defined by effect/NES direction, not by significance overlap. Subset-only FDR hits are highlighted separately.</div>']
    if not summary.empty:
        H+=['<h2>Cell-type overview</h2>',html_table(summary,list(summary.columns),100)]
    gene_cols=['gene',f'{a.reference_analysis}_log2FoldChange',f'{a.reference_analysis}_FDR']+[f'{x}_log2FoldChange' for x in a.analyses if x!=a.reference_analysis]+['n_significant','direction_class','median_abs_delta_from_reference']
    path_cols=['Term','gene_set',f'{a.reference_analysis}_NES',f'{a.reference_analysis}_FDR q-val']+[f'{x}_NES' for x in a.analyses if x!=a.reference_analysis]+['n_significant','direction_class','median_abs_delta_from_reference']
    for cell in cells:
        H.append(f'<h2>{html.escape(cell)}</h2>')
        g=genes[genes.cell_type==cell] if not genes.empty else pd.DataFrame(); p=paths[paths.cell_type==cell] if not paths.empty else pd.DataFrame()
        if not g.empty:
            H.append('<h3>Concordant reference DE hits</h3>'); H.append(html_table(g[g.standout_type=='REFERENCE_HIT_CONCORDANT'].sort_values('agreement_priority',ascending=False),gene_cols,a.top_n))
            H.append('<h3>Genes requiring attention</h3>'); H.append(html_table(g[g.standout_type.isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT','SUBSET_SPECIFIC_HIT'])].sort_values(['n_significant','mean_abs_log2FC'],ascending=False),gene_cols,a.top_n))
        if not p.empty:
            H.append('<h3>Concordant reference pathways</h3>'); H.append(html_table(p[p.standout_type=='REFERENCE_HIT_CONCORDANT'].sort_values('agreement_priority',ascending=False),path_cols,a.top_n))
            H.append('<h3>Pathways requiring attention</h3>'); H.append(html_table(p[p.standout_type.isin(['REFERENCE_HIT_DISCORDANT','DIRECTION_DISCORDANT','SUBSET_SPECIFIC_HIT'])].sort_values(['n_significant','mean_abs_NES'],ascending=False),path_cols,a.top_n))
    H.append('</body></html>'); (out/'interpretation_report.html').write_text('\n'.join(H),encoding='utf-8')
    print(f'Reports written to: {out}')

if __name__=='__main__': main()
