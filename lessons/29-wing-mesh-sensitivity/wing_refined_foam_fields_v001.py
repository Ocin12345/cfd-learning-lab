"""Read-only native checks for the saved 3D finite-wing Foundation14 run.

--run names the actual OpenFOAM case, --out a fresh external analysis folder.
Exit2 retains genuine missed numerical targets. No CFD or WSL is started.
"""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
import numpy as np
from prepare_openfoam_finite_wing_v003 import (
    PATCHES, FLOAT, require, sha, clean, dictionaries, native_poly, native_patches,
    mesh_quality, json_new, copied,
)


class Evidence:
    def __init__(self):
        self.before = {}

    def text(self, path):
        path = Path(path).resolve()
        require(path.is_file(), 'Missing actual native evidence: ' + str(path))
        if path not in self.before:
            self.before[path] = sha(path)
        return path.read_text(encoding='utf-8-sig')

    def read_json(self, path):
        return json.loads(self.text(path))

    def finish(self):
        rows = [dict(path=str(path),bytes=path.stat().st_size,sha256_before=before,sha256_after=sha(path))
                for path,before in self.before.items()]
        require(all(r['sha256_before']==r['sha256_after'] for r in rows), 'Native source changed during read-only analysis')
        return dict(inputs_unchanged_during_analysis=True,sources=rows,solver_started_or_modified=False)


def values(raw, count, components):
    text = clean(raw)
    require(re.search(r'\bformat\s+ascii\s*;',text), 'Native ASCII field required')
    require(not re.search(r'(?i)(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])',text), 'Nonfinite native field/boundary token')
    m = re.search(r'\binternalField\s+uniform\s+([^;]+);',text)
    if m:
        row = np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
        require(row.size==components and np.isfinite(row).all(),'Invalid uniform native field')
        return np.tile(row,(count,1)) if components>1 else np.full(count,row[0])
    kind = 'vector' if components==3 else 'scalar'
    m = re.search(r'\binternalField\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count, 'Native internal-field count differs')
    data = np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ')
    require(data.size==count*components and np.isfinite(data).all(),'Native field length/finite check failed')
    return data.reshape((count,components)) if components>1 else data


def field(raw,count,components,dimensions):
    d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',clean(raw))
    require(d is not None and list(map(int,d[1].split()))==dimensions,'Native field dimensions differ')
    return values(raw,count,components)


def boundary_mass(raw,patches):
    text=clean(raw)
    require('boundaryField' in text,'Missing phi boundary fields')
    blocks=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',text.split('boundaryField',1)[1]))
    totals={}
    inward=0.0
    outward=0.0
    for name in PATCHES:
        require(name in blocks,'Missing native phi patch '+name)
        body=blocks[name]
        m=re.search(r'\bvalue\s+uniform\s+('+FLOAT+r')\s*;',body)
        if m:
            data=np.full(patches[name]['count'],float(m[1]))
        else:
            m=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',body,re.S)
            require(m is not None and int(m[1])==patches[name]['count'],'Native phi patch count differs '+name)
            data=np.fromstring(m[2],sep=' ')
        require(data.size==patches[name]['count'] and np.isfinite(data).all(),'Nonfinite/wrong phi boundary length')
        totals[name]=float(data.sum())
        inward-=float(data[data<0].sum())
        outward+=float(data[data>0].sum())
    require(inward>0,'Missing native inward boundary flux')
    net=sum(totals.values())
    return dict(signed_outward_volume_flux_m3_s=totals,total_inward_volume_flux_m3_s=inward,
        total_outward_volume_flux_m3_s=outward,net_outward_volume_flux_m3_s=net,
        relative_imbalance=abs(net)/inward,
        definition='abs(sum outward phi over every physical boundary face)/sum magnitude of inward face phi; includes farfield and wing')


def residual_history(text):
    require(not re.search(r'FOAM FATAL|Segmentation fault|Floating point exception\s*$',text,re.M),'Actual solver fatal error')
    sections=re.split(r'^Time = ('+FLOAT+r')s?\s*$',text,flags=re.M)
    rows=[]
    for i in range(1,len(sections),2):
        iteration=float(sections[i])
        samples={name:[] for name in ('p','Ux','Uy','Uz')}
        for name,initial in re.findall(r'Solving for (p|Ux|Uy|Uz), Initial residual = ('+FLOAT+r')',sections[i+1]):
            samples[name].append(float(initial))
        require(all(samples.values()),'Missing actual p/Ux/Uy/Uz residual at iteration '+str(iteration))
        row=dict(iteration=int(iteration),**{name:max(v) for name,v in samples.items()})
        require(iteration==row['iteration'] and all(np.isfinite(v) for v in row.values()),'Invalid steady residual samples')
        rows.append(row)
    require([r['iteration'] for r in rows]==list(range(1,1501)),'Need actual fixed-cap iterations1..1500')
    return rows


def forces(text):
    header=re.search(r'^#\s*Time\s+Cm\s+Cd\s+Cl\s+Cl\(f\)\s+Cl\(r\)\s*$',text,re.M)
    require(header is not None,'Native forceCoeffs columns differ from Foundation14 official source')
    rows=[]
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        data=list(map(float,line.split()))
        require(len(data)==6 and all(np.isfinite(data)),'Invalid native force row')
        iteration=data[0]
        require(iteration.is_integer(),'Steady force index must be an iteration, not physical time')
        if iteration>0:
            rows.append(dict(iteration=int(iteration),Cm=data[1],Cd=data[2],Cl=data[3]))
    indices=[r['iteration'] for r in rows]
    require(indices==sorted(set(indices)) and all(1<=i<=1500 for i in indices),
            'Native force iteration indices must be unique, increasing, within1..1500')
    return rows

