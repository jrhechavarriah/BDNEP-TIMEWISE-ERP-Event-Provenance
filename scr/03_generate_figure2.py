#!/usr/bin/env python3
"""Rebuild main Figure 2 only from hash-verified P14 waveform CSV, no frozen contrasts as input."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from bdnep_inputs import verified_artifact
import numpy as np

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verified(replay, name):
 return verified_artifact(replay, name)

def main():
 a=argparse.ArgumentParser();inputs=a.add_mutually_exclusive_group(required=True);inputs.add_argument('--repo',dest='replay_dir',type=Path);inputs.add_argument('--replay-dir',dest='replay_dir',type=Path);a.add_argument('--out-dir',type=Path,required=True);p=a.parse_args()
 report=json.loads(verified(p.replay_dir,'replay_verification_report.json').read_text())
 if report.get('status')!='FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES' or report.get('n_subjects')!=37:raise RuntimeError('Replay not verified')
 with verified(p.replay_dir,'figure_01_PO8_waveform_data.csv').open(newline='',encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
 t=np.array([float(x['time_ms']) for x in rows]);o=np.array([float(x['original_face_minus_car_uV']) for x in rows]);c=np.array([float(x['corrected_face_minus_car_uV']) for x in rows]);
 if len(t)<900 or not np.all(np.diff(t)>0) or not all(np.isfinite(v).all() for v in (t,o,c)):raise RuntimeError('Bad waveform values')
 # Frozen ROI uses actual discrete grid; 41 samples between 110 and 150 inclusive.
 roi=(t>=110-1e-9)&(t<=150+1e-9)
 if roi.sum()!=41 or abs((c[roi]-o[roi]).mean()-report['group']['PO8_A_mean_delta_uV'])>1e-6:raise RuntimeError('ROI differs from validated group summary')
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 if p.out_dir.exists():raise RuntimeError('Refusing to overwrite existing figure directory')
 p.out_dir.mkdir(parents=True)
 fig,ax=plt.subplots(figsize=(9,5))
 ax.plot(t,o,label='Original digital markers');ax.plot(t,c,label='Author-corrected digital markers');ax.axvspan(110,150,alpha=.14,label='Prespecified ROI');ax.axhline(0,color='gray',lw=.7)
 ax.set(xlabel='Time relative to digital marker (ms)',ylabel='Face minus car, PO8 (µV)',title='Grand-average PO8 ERP contrast (n = 37)')
 ax.set_xlim(-200,800);ax.legend(fontsize=8)
 fig.tight_layout();dst=p.out_dir/'Figure_2_PO8_ERP.png';fig.savefig(dst,dpi=300);fig.savefig(dst.with_suffix('.svg'));plt.close(fig)
 print('FIGURE_2_GENERATED',dst);print('SHA256',digest(dst))
if __name__=='__main__':main()
