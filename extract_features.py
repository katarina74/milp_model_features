from multiprocessing import Process
import ctypes
import pandas as pd
from os import listdir
from os.path import isfile, join
import re
from tqdm import tqdm
from os import listdir
from os.path import isfile, join


def write_fetures_pre(lp_file, collection_dir, collection_csv_dir):
    instance_name = re.split(r'[.]', lp_file)[0]
    
    mf.write_features(bytes(f'{collection_dir}/{lp_file}', encoding='utf8'), 
                      bytes(f'{collection_csv_dir}/{instance_name}.csv', encoding='utf8'),
                      bytes(instance_name, encoding='utf8'))

mf = ctypes.CDLL(".ModelFeatures/x64/Debug/ModelFeatures.dll")
lp_files = [f for f in listdir("collection") if isfile(join("collection", f))]

if __name__ == '__main__':
    for lp_file in tqdm(lp_files):
        p1 = Process(target=write_fetures_pre, args=(lp_file, "collection", "collection_csv"))
        p1.start()
        p1.join(10*60)
        if p1.is_alive(): p1.terminate()


    csv_files = [f for f in listdir("collection_csv") if isfile(join("collection_csv", f))]

    df_list = []
    for csv_file in csv_files:
        df = pd.read_csv(f"collection_csv/{csv_file}").set_index('INSTANCE_NAME').T
        df_list.append(df)
    df_from_scip = pd.concat(df_list)

    df_from_miplib = pd.read_csv("TheCollectionSet.csv", index_col="Instance  Ins.")
    data = pd.merge(df_from_scip, df_from_miplib, left_index=True, right_index=True)

    data.to_csv("miplib.csv")
    