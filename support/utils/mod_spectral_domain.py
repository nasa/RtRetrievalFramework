#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# modifies spectral domain (through selection or masking) specified in hdf5 input file
# by F. Oyafuso

from __future__ import print_function
from builtins import str
from builtins import range
import os, sys, shutil, re
import numpy, scipy, h5py
import shlex, subprocess
import argparse

parser = argparse.ArgumentParser(description='parser for mod_spectral_domain.py', prefix_chars='-')
parser.add_argument('--band', dest="band", nargs=1, default=['WCO2'],
                    help='∊ {WCO2,SCO2,O2A}')
parser.add_argument('--lines', dest="line_method", nargs=1, default=['mask'],
                    help='∊ {mask,select}')
parser.add_argument('--linetype', dest="line_type", nargs=1, default=['solar'],
                    help='∊ {solar,earth}')
parser.add_argument('--halfwidth', dest="win_halfwidth_dflt", nargs=1, default=[0.2],
                    help='default value of window halfwidth')
parser.add_argument('--linestrength', dest="strength_cutoff", nargs=1, default=[7e-24],
                    help='cutoff of linestrength')
parser.add_argument('--solar_od', dest="solar_od_cutoff", nargs=1, default=[0.1],
                    help='cutoff of solar optical depth')
parser.add_argument('--in', dest="fil_in", nargs=1, default=[os.environ['L2_INPUT_PATH']+'/fts/input/l2_fts_static_input.h5'],
                    help='input file')
parser.add_argument('--out', dest="fil_out", nargs=1, default=['./l2_fts_static_input.h5'],
                    help='output file')
parser.add_argument('--linelist', dest="fil_ll", nargs=1, default=[os.environ['L2_ABSCO_PATH']+'/input/line_lists/HITRAN08.par'],
                    help='output file')
args = parser.parse_args()

window_dflt={}
window_dflt['O2A']=[12950,13200]
window_dflt['WCO2']=[6173,6275]
window_dflt['SCO2']=[4785,4925]
band_dflt={}
band_dflt['O2A']=0
band_dflt['WCO2']=1
band_dflt['SCO2']=2

w=window_dflt[ args.band[0] ]
band=band_dflt[ args.band[0] ]

if args.line_type[0]=='solar':
    print('Reading \'' + args.fil_in[0] + '\'')
    fp1=h5py.File(args.fil_in[0],'r')
    frqs0=scipy.array(fp1['/Solar/Solar_Line_List/frequency'])
    ss0=scipy.array(fp1['/Solar/Solar_Line_List/optical_thickness'])
    fp1.close()

    idx=(ss0>args.solar_od_cutoff[0]) & (frqs0>w[0]) & (frqs0<w[1])
    frqs=frqs0[idx]
    win_halfwidth=scipy.ones(frqs.size) * args.win_halfwidth_dflt[0]
elif args.line_type[0]=='earth':
    print('Reading \'' + args.fil_ll[0] + '\'')
    cmd=['perl','-lne',
         ('$nu=substr($_,3,12); $S=substr($_,16,9); print "$nu" if /^ 2/ && {0}<$nu && $nu<{1} && $S>'+str(args.strength_cutoff[0])+';').format(w[0],w[1]),
         args.fil_ll[0]]
    #dum=subprocess.check_output(cmd)
    p=subprocess.Popen(cmd, stdout=subprocess.PIPE)
    frqs_str, err = p.communicate()
    frqs=[float(s) for s in frqs_str.split()]
    win_halfwidth=scipy.ones(frqs.size) * args.win_halfwidth_dflt[0]

# block out additional line to block out
if args.band[0]=='WCO2':
    line=6250.4
    idx=scipy.where( (line < frqs) & (frqs < w[1]))
    if idx[0].size>0:
        frqs=scipy.insert(frqs, idx[0][0], line)
        win_halfwidth=scipy.insert(win_halfwidth, idx[0][0], 0.15)

print('Copying \''+args.fil_in[0]+'\' to \''+args.fil_out[0]+'\'')
shutil.copy2(args.fil_in[0],args.fil_out[0])

fp=h5py.File(args.fil_out[0],'a')

mw_orig=scipy.array(fp['/Spectral_Window/microwindow'])

Nspec=3
if args.line_method[0]=='select':
    Nwin=len(frqs)
    mw_test=scipy.zeros([Nspec,Nwin,2], dtype=scipy.float64)
    # leftmost window
    mw_test[band,0,0]=frqs[0] - win_halfwidth[0]
    if (mw_test[band,0,0] < w[0]):
        mw_test[band,0,0]=w[0]
    mw_test[band,0,1]=frqs[0] + win_halfwidth[0]
    # middle windows
    counter=1
    for j in range(1,Nwin):
        edge_left=frqs[j] - win_halfwidth[j]
        if edge_left > mw_test[band,counter-1,1]:
            mw_test[band,counter,0]=edge_left
            mw_test[band,counter,1]=frqs[j] + win_halfwidth[j]
            counter=counter+1
        else:
            # expand previous window instead of adding another one
            print('warning: expanding pervious window: '+str(mw_test[band,counter-1,1])+' -> '+str(frqs[j] + win_halfwidth[j]))
            mw_test[band,counter-1,1]=frqs[j] + win_halfwidth[j]
    # check rightmost window
    if (mw_test[band,counter-1,1] > w[1]):
        mw_test[band,counter-1,1]=w[1]
    mw_new=mw_test[:,0:counter,:]
elif args.line_method[0]=='mask':
    Nwin=len(frqs)+1
    mw_test=scipy.zeros([Nspec,Nwin,2], dtype=scipy.float64)
    mw_test[band,0,0]=mw_orig[band,0,0]
    # middle windows
    counter=1
    for j in range(1,Nwin):
        edge_right=frqs[j-1] - win_halfwidth[j-1]
        if edge_right > mw_test[band,counter-1,0]:
            mw_test[band,counter-1,1]=edge_right
            mw_test[band,counter,0]=frqs[j-1] + win_halfwidth[j-1]
            counter=counter+1
        else:
            print('warning: deleting window starting at '+str(mw_test[band,counter-1,0]))
            mw_test[band,counter-1,0]=frqs[j-1] + win_halfwidth[j-1]
    # rightmost window
    mw_test[band,counter-1,1]=mw_orig[band,0,1]
    mw_new=mw_test[:,0:counter,:]

print([counter, Nwin])
print(frqs)
print(win_halfwidth)
print(mw_new[band,:,:])

attr=fp['/Spectral_Window/microwindow'].attrs['Units']

del fp['/Spectral_Window']['microwindow']

fp['/Spectral_Window/microwindow']=mw_new
fp['/Spectral_Window/microwindow'].attrs['Units']=attr

if False:
    ap0=scipy.array(fp['/Instrument/Continuum/a_priori'])
    ap0[2,0]=1.97526
    ap0[2,1]=-3.44291e-4
    del fp['/Instrument/Continuum']['a_priori']
    fp['/Instrument/Continuum/a_priori']=ap0

fp.close()

