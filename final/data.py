import os, glob
import pandas as pd
import numpy as np

def readPSSMFile(inputFilePath, SeqName=False):    
    startCol = 2; 
    if SeqName == True:
        startCol = 1; 
    # 22 amino acid matrix, 42 means with their frequence
    endCol = 22; 
    # Set some of unuseful rows by number of row 
    setSkipRows = (0,1,2);
    # Read data
    pssmData = pd.read_csv(inputFilePath, 
                            skiprows=setSkipRows,
                            # White space delimater ignores more spaces
                            delim_whitespace=True,
                            # No header
                            header=None,
                            usecols=range(startCol, endCol)
                            )
    # Remove last 5 rows
    maxRow = (len(pssmData.index)) - 6 
    # Select only content of PSSM profiles
    pssmData = pssmData.loc[0:maxRow, :]

    return pssmData

def make_pssm_dir(pssm_path):
    # Make directory to store
    # os.makedirs("train_pssm", exist_ok=True)
    dir_path = pssm_path.split('.')[0] + '_' + pssm_path.split('.')[-1]
    os.makedirs(dir_path, exist_ok=True)
    with open(pssm_path) as f:
        data = f.readlines()
    str_idx_list = [] # [,]
    idx_list = [] # Query profile of sequence *
    for ii, i in enumerate(data):
        if i.split(' ')[0] == "Query":
            idx = i.split(' ')[-1][:-1]
            idx_list.append(idx)
            str_idx_list.append(ii)
    # print((str_idx_list))
    for i,j in zip(range(len(str_idx_list)), idx_list):
        p = os.path.join(dir_path,str(j))+".pssm"
        if i != (len(str_idx_list)-1): 
            with open(p, 'w') as f:
                f.write(''.join(data[str_idx_list[i]:str_idx_list[i+1]]))
        elif i == (len(str_idx_list)-1):
            with open(p, 'w') as f:
                f.write(''.join(data[str_idx_list[i]:len(data)]))


def extract_data(path):
    # Let pssm matrix  normalize to numpy 
    dir_path = path.split('_')[0] + '_' + "np"
    # path = "train_pssm_test"
    os.makedirs(dir_path, exist_ok=True)
    am = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    pssm_paths = glob.glob(os.path.join(path, '*'))
    # print(pssm_paths)
    for pssm_path in pssm_paths:
        df = readPSSMFile(pssm_path, True)
        # print(df)
        # final_array = np.empty((20,20))
        final_array = []
        for idx,i in enumerate(am):
            # print(i)
            # print(df[1])
            s_df = df[df[1] == i]
            if len(s_df) == 0:
                sigmoid_array = np.zeros(20)
            else:
                # print("??????????????", len(s_df))
                # print(":::::;",s_df)
                # s_df = df.iloc[:,1:]
                s_df = s_df.iloc[:,1:]
                # if idx == len(am)-1:
                    # print(s_df)
                # print(s_df)
                s_df_mean = s_df.mean()
                np_data = s_df_mean.to_numpy()
                # print(np_data)
                # print(np_data.shape)
                sigmoid_array = 1 / (1 + np.exp(-np_data))
            final_array.append(sigmoid_array)
            # print(sigmoid_array)
        final_array = np.array(final_array)
        # print(final_array)
        with open(os.path.join(dir_path, pssm_path.split('/')[-1].split('.')[0]+'.npy'), 'wb') as f:
            np.save(f, final_array)
        # print(final_array.shape)

def main():
    make_pssm_dir("train.pssm")
    make_pssm_dir("test.pssm")
    extract_data("train_pssm")
    extract_data("test_pssm")


if __name__ == '__main__':
    main()
