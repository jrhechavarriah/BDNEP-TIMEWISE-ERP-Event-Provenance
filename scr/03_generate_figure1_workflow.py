#!/usr/bin/env python3
"""Create a SOURCE-VERIFIED editorial schematic of digital event-marker provenance.

Not a neural measurement or signal estimate. The plot is explanatory rather
than a recomputation of EEG waveforms; physical screen onset is not measured.
The original editorial v3 may use a different design and needs separate review.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

EXPECTED_CONTRACT_SHA='f75505d2b30f57d4b6ce88a5efcd60c2cf27398e7021dbe2545f34dab4d4827c'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--contract',type=Path,required=True)
 p.add_argument('--replay-report',type=Path,required=True)
 p.add_argument('--out-dir',type=Path,required=True)
 a=p.parse_args()
 if sha(a.contract)!=EXPECTED_CONTRACT_SHA:raise RuntimeError('Source manifest has changed')
 c=json.loads(a.contract.read_text(encoding='utf-8-sig'))
 r=json.loads(a.replay_report.read_text(encoding='utf-8-sig'))
 if c.get('file_count')!=261 or len(c.get('cohort',[]))!=37:raise RuntimeError('Source cohort mismatch')
 if r.get('status')!='FULL_INDEPENDENT_REPLAY_MATCHES_FROZEN_NUMERICAL_REFERENCES' or r.get('n_subjects')!=37:raise RuntimeError('Replay not verified')
 if r.get('sample_event_shift')!={'stimuli':26,'responses':0} or r.get('orphan_responses')!=67:raise RuntimeError('Digital event lock mismatch')
 if r['group']['n_A_stimuli']!=11840 or r['group']['n_B_stimuli']!=10810:raise RuntimeError('Group totals mismatch')
 if a.out_dir.exists():raise RuntimeError('Refusing to overwrite figure output')
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import FancyBboxPatch
 a.out_dir.mkdir(parents=True)
 fig,ax=plt.subplots(figsize=(12.2,6.4));ax.set_xlim(0,12.2);ax.set_ylim(0,6.4);ax.axis('off')
 def box(x,y,w,h,title,lines):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.10',lw=1.3,fill=False))
  ax.text(x+w/2,y+h-.40,title,ha='center',va='center',weight='bold',fontsize=12)
  ax.text(x+w/2,y+h*.44,lines,ha='center',va='center',fontsize=9.6,linespacing=1.5)
 box(.3,3.35,3.0,2.1,'OSF event records',
     'Original stimulus markers\n11,840 stimuli + 11,814 responses\n37 participants')
 box(4.55,3.35,3.0,2.1,'Author-corrected records',
     'Stimulus identities preserved\nVisual stimuli: +26 samples\nResponses: unchanged')
 box(8.75,3.35,3.0,2.1,'NEMAR N170 EEG',
     'The same continuous recordings\n1024 Hz; 30-channel EEG reference\nNo artifact rejection')
 box(2.05,.45,8.1,1.7,'Frozen paired ERP analysis',
     'Identical signal filtering and epoching under both marker references\nPO8 face minus car, 110–150 ms (primary)\nPO7 and behaviorally restricted Arm B (descriptive)')
 ax.annotate('',xy=(4.45,4.35),xytext=(3.40,4.35),arrowprops={'arrowstyle':'-|>','lw':1.7})
 ax.annotate('',xy=(8.65,4.35),xytext=(7.65,4.35),arrowprops={'arrowstyle':'-|>','lw':1.7})
 ax.plot([10.25,10.25,6.1],[3.25,2.76,2.76],lw=1.6)
 ax.annotate('',xy=(6.1,2.25),xytext=(6.1,2.76),arrowprops={'arrowstyle':'-|>','lw':1.7})
 ax.text(6.1,5.85,'Digital marker comparison: +26 samples = 25.390625 ms at 1024 Hz',ha='center',fontsize=12,weight='bold')
 ax.text(6.1,.14,'Schematic only — physical stimulus onset, display latency and timing jitter were not measured.',ha='center',fontsize=9,style='italic')
 fig.tight_layout()
 out=a.out_dir/'Figure_1_Digital_Event_Marker_Workflow.png'
 fig.savefig(out,dpi=300,bbox_inches='tight')
 fig.savefig(out.with_suffix('.svg'),bbox_inches='tight')
 plt.close(fig)
 print('FIGURE1_EDITORIAL_CANDIDATE_GENERATED:',out)
 print('NO_ONSET_CALIBRATION_PERFORMED')
if __name__=='__main__':main()
