import pandas as pd
from io import StringIO

po_po_table = """
AA,-0.83,-8.10,-22.5
UU,-0.83,-8.10,-22.5
AU,-0.56,-5.53,-15.4
UA,-0.58,-6.40,-18.0
CA,-0.95,-6.89,-18.4
UG,-0.95,-6.89,-18.4
GU,-0.94,-7.12,-19.1
AC,-0.94,-7.12,-19.1
CU,-0.94,-7.51,-20.3
AG,-0.94,-7.51,-20.3
GA,-0.88,-6.51,-17.4
UC,-0.88,-6.51,-17.4
CG,-1.62,-10.81,-28.5
GC,-1.76,-12.68,-33.8
GG,-1.09,-6.09,-15.5
CC,-1.09,-6.09,-15.5
EA,0.49,20.73,62.6
UE,0.49,20.73,62.6
AE,0.48,20.20,61.0
EU,0.48,20.20,61.0
EC,0.63,21.07,63.2
GE,0.63,21.07,63.2
CE,0.43,18.09,54.7
EG,0.43,18.09,54.7
"""



df_po_po = pd.read_csv(
    StringIO(po_po_table),
    header=None,
    names=['nucleotide', 'G_50', 'H', 'S'],
)

df_po_po['G_37'] = (df_po_po['H'] * 1000 - 310.15 * df_po_po['S']) / 1000
po_po_g37_raw = df_po_po.set_index('nucleotide')['G_37']
po_po_g37 = po_po_g37_raw.to_dict()
po_po_g50_raw = df_po_po.set_index('nucleotide')['G_50']
po_po_g50 = po_po_g50_raw.to_dict()


ps_po_table = """
AA,-0.52,-5.81,-16.4
AU,-0.42,-5.64,-16.2
AC,-0.88,-8.66,-24.1
AG,-0.71,-6.13,-16.8
UA,-0.30,-3.86,-11.0
UU,-0.49,-5.87,-16.7
UC,-0.64,-6.13,-17.0
UG,-0.77,-7.28,-20.2
CA,-0.82,-7.23,-19.8
CU,-0.66,-6.30,-17.5
CC,-0.91,-5.57,-14.4
CG,-1.21,-8.07,-21.2
GA,-0.85,-7.92,-21.9
GU,-0.61,-5.46,-15.0
GC,-1.15,-6.63,-17.0
GG,-1.09,-7.75,-20.6
EA,-0.57,13.64,44.0
AE,-0.54,15.07,48.3
EU,-0.59,14.87,47.8
UE,-0.56,14.72,47.3
EC,-0.58,10.31,33.7
CE,-0.55,10.49,34.2
EG,-0.35,16.13,51.0
GE,-0.44,14.65,46.7
"""



df_ps_po = pd.read_csv(
    StringIO(ps_po_table),
    header=None,
    names=['nucleotide', 'G_50', 'H', 'S'],
)


df_ps_po['G_37'] = (df_ps_po['H'] * 1000 - 310.15 * df_ps_po['S']) / 1000
ps_po_g37_raw = df_ps_po.set_index('nucleotide')['G_37']
ps_po_g37 = ps_po_g37_raw.to_dict()
ps_po_g50_raw = df_ps_po.set_index('nucleotide')['G_50']
ps_po_g50 = ps_po_g50_raw.to_dict()


ps_diff_g37 = (ps_po_g37_raw - po_po_g37_raw).to_dict()
ps_diff_g50 = (ps_po_g50_raw - po_po_g50_raw).to_dict()

for key,value in ps_diff_g50.items():
    print(f"{key},   {100 * value}")