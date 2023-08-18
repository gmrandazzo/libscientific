import sys
import seaborn as sns
import csv
import pandas as pd

def read_csv_file(fcsv):
    return pd.read_csv(fcsv, header=0)


def plot(d, savefname="input_vs_something.png", y_label='CPU time (sec)'):
    p =sns.lineplot(data=d)
    p.get_legend().set_visible(False)
    p.set_xlabel('input size (instances x features)', fontdict={'size': 15})
    p.set_ylabel(y_label, fontdict={'size': 15})
    p.get_figure().savefig(savefname)
    
def prepare_input_size_vs_cpu_time(ds):
    tp = {}
    indx = []
    h = []
    for key in ds.keys():
        if len(tp) == 0:
            key_ = key.replace(".csv", "")
            h.append(key_)
            tp[key_] = []
            for i in range(len(ds[key])):
                indx.append(ds[key]["row"][i]*ds[key]["cols"][i])
                tp[key_].append(ds[key]["time(sec)"][i])
        else:
            key_ = key.replace(".csv", "")
            h.append(key_)
            tp[key_] = []
            for i in range(len(ds[key])):
                tp[key_].append(ds[key]["time(sec)"][i])
    return pd.DataFrame(tp, index=indx, columns=h)

ds = {}
for i in range(1, len(sys.argv)):
    ds[sys.argv[i]] = read_csv_file(sys.argv[i])

plot(prepare_input_size_vs_cpu_time(ds), "input_vs_cputime.png", 'CPU time (sec)')

    
            

